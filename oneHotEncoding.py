import sys, numpy as np, random
import formula_processing as fp

def argparsr(argl):
    ##argl is the list of arguments the user inputted
    for ar in range(1, len(argl)):
        print(ar)
        if argl[ar] == '-inp':
            inpf = argl[ar + 1]
        elif argl[ar] == '-binw':
            binw = argl[ar + 1]
        elif argl[ar] == '-output':
            outp = argl[ar + 1]
        elif argl[ar] == '-perct':
            perc = argl[ar + 1]
        elif argl[ar] == '-makeR':
            makr = argl[ar + 1]
        elif argl[ar] == '-makeA':
            maka = argl[ar + 1]
        elif argl[ar] == '-forms':
            frm = argl[ar + 1]
    return(inpf, binw, outp, perc, makr, maka, frm)

def show_help():
    print('ARGS:\n'\
        '-inp full path to input file\n'\
        '-output full path to desired output file\n'\
        '-perct what percentage (as decimal) do you want for the test set?\n'\
        '-makeR do you want to make iRF input? y OR n\n'\
        '-makeA do you want to make WEKA input? y OR n\n'\
        '-forms do you want to include formula features as well?\n'\
        '-binw what is your preferred binwidth?')

def parse_mgf(inp, frm):

    ##inp is a file object of the input mgf file
    ##Frm is a string (y OR n)
    ##RETURNS: Dict, with class as key and name,
    ##m/z, and MS/MS peaks as values.
    mgfd = {}
    strt = True
    for inpl in inp:
        if inpl.startswith('Name') and strt == True:
            nam = inpl[6:].strip()
            clas = nam.split()[0].replace('[','').replace(']','')
            name = nam.replace(' ', '_').replace(':', '-').replace(';', '').replace('}','').replace(',','-')
            mgfd[clas] = []
            strt = False

        elif inpl.startswith('Name') and strt != True:
            if clas in mgfd.keys():
                if frm == 'y':
                    mgfd[clas].append([name, clas, mz, peaks, form])
                else:
                    mgfd[clas].append([name, clas, mz, peaks])
            else:
                if frm == 'y':
                    mgfd[clas] = [[name, clas, mz, peaks]]
                else:
                    mgfd[clas] = [[name, clas, mz, peaks]]

            ##resetting variables for next entry
            nam = inpl[6:].strip()
            clas = nam.split()[0].replace('[','').replace(']','')
            name = nam.replace(' ', '_').replace(':', '-').replace(';','').replace('}','').replace(',', '-')

        elif inpl.startswith('MW:'):
            mz = float(inpl.strip().split(': ')[1])
        elif inpl.startswith('Comment'):
            form = inpl.strip().split(';')[-2].strip()
        elif inpl.startswith('Num Peaks:'):
            peaks = []
            inpl =inp.readline()
            while inpl[:1].isdigit():
                ##we are just logging mz, not abund
                peaks.append(float(inpl.strip().split()[0]))
                inpl = inp.readline()

    mgfd[clas].append([name, clas, mz, peaks])
    return(mgfd)

def onehot(lin, bigl):
    ##lin is a list with a sublist of MS/MS peaks
    ##bigl is the list of ranges for encoding
    ##RETURNS: lin with the sublist of peaks
    ##one-hot encoded.
    peaks = lin[3]
    ohpeaks = []
    for ran in range(len(bigl)-1):
        st = float(bigl[ran])
        end = float(bigl[ran + 1])
        if (any(peak >= st and peak <= end for peak in peaks)):
            ohpeaks.append(1)
        else: ohpeaks.append(0)

    lin[3] = ohpeaks
    return lin

def make_arff(bnw, mgfd, out, p, frm):
    ##mgfd is the dict containing output of parse_mgf
    ##binwidth is float or int
    ##OUTPUTS: arff file with peaks as
    ##one-hot encoded features

    #Initialize ARFF file
    train=open(out + 'training.arff','w')
    test=open(out + 'test.arff','w')

    biglist=np.arange(50,2500,float(bnw))
    train.write('@relation\tMSMS_multiclass\n\n')
    test.write('@relation\tMSMS_multiclass\n\n')

    train.write('@attribute\tID\tstring\n')
    test.write('@attribute\tID\tstring\n')
    
    for i in range(len(biglist)-1):
        train.write(f'@attribute\t{biglist[i]}-{biglist[i+1]}\tnumeric\n')
        test.write(f'@attribute\t{biglist[i]}-{biglist[i+1]}\tnumeric\n')

    if frm == 'y':
        frmft = fp.make_arff_header()
        test.write(f'{frmft}')
        train.write(f'{frmft}')

    ##adding line for all possible classes
    clssls = mgfd.keys()

    clsstr = ','.join(clssls)
    print(clsstr)
    train.write(f'@attribute\tclass\t{{{clsstr}}}\n\n')
    test.write(f'@attribute\tclass\t{{{clsstr}}}\n\n')
    train.write('@data\n')
    test.write('@data\n')

    ##reserving instances for the test arff
    rlst = []
    trlst = []
    for insl in mgfd.values():
        if len(insl) >= 10:
            sampnum = int(float(p) * len(insl))
        else:
            sampnum = 1
        tmp = random.sample(insl, sampnum)

        ##now getting remainder
        remainder = insl
        for x in tmp:
            remainder.remove(x)

        ##getting instances not in tmp (these are the training ones)
        rlst.append(tmp)
        trlst.append(remainder)
        
    flatr = [item for sublist in rlst for item in sublist]
    trainr = [item for sublist in trlst for item in sublist]
    print('done with making lists!')

    ##doing one-hot now and writing out
    ##first for test arff instances
    #for ind in range(len(flatr)):
    for ind in range(3):
        flatr[ind] = onehot(flatr[ind], biglist)
        featstr = ','.join(str(x) for x in flatr[ind][3])
        if frm == 'y':
            c, h, o, n, p, s, mim = fp.parse_form(flatr[ind][4])
            h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
            mad, amd, rmd = fp.mass_defects(mim)
            ffeatstr = ','.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
             h2o, c2p, c2n, c2s, amd, mad, rmd])
            test.write(f'{flatr[ind][0]},{featstr},{ffeatstr},{flatr[ind][1]}\n')
        else:
            test.write(f'{flatr[ind][0]},{featstr},{flatr[ind][1]}\n')
    test.close()
    print('Done with test!')

    ##now training
    #for ind in range(len(trainr)):
    for ind in range(3):
        trainr[ind] = onehot(trainr[ind], biglist)
        featstr = ','.join(str(x) for x in trainr[ind][3])
        if frm == 'y':
            c, h, o, n, p, s, mim = fp.parse_form(trainr[ind][4])
            h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
            mad, amd, rmd = fp.mass_defects(mim)
            ffeatstr = ','.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
             h2o, c2p, c2n, c2s, amd, mad, rmd])
            train.write(f'{trainr[ind][0]},{featstr},{ffeatstr},{trainr[ind][1]}\n')
        else: 
            train.write(f'{trainr[ind][0]},{featstr},{trainr[ind][1]}\n')
    train.close()
    print('Done with train!')

def make_irf(bnw, mgfd, out, p, frm):
    ##mgfd is the dict containing output of parse_mgf
    ##binwidth is float or int
    ##OUTPUTS: x: one hot features for each instance, 
    # y: class of each instance, xtest: features for
    # each test instance, ytest: class of each test ins

    trainx = open(out + 'xtrain.tab', 'w')
    testx = open(out + 'xtest.tab', 'w')
    trainy = open(out + 'ytrain.tab', 'w')
    testy = open(out + 'ytest.tab', 'w')

    biglist=np.arange(50,2500,float(bnw))

    ##writing out feature names for x-files
    news = ''
    for i in range(len(biglist)-1):
        news += f'{biglist[i]}-{biglist[i+1]}\t'
    
    if frm == 'y':
        news += fp.make_irf_header()

    trainx.write(f'{news}\n')
    testx.write(f'{news}\n')

    ##reserving instances for the test arff
    rlst = []
    trlst = []
    for insl in mgfd.values():
        if len(insl) >= 10:
            sampnum = int(float(p) * len(insl))
        else:
            sampnum = 1
        tmp = random.sample(insl, sampnum)

        ##now getting remainder
        remainder = insl
        for x in tmp:
            remainder.remove(x)

        ##getting instances not in tmp (these are the training ones)
        rlst.append(tmp)
        trlst.append(remainder)
        
    flatr = [item for sublist in rlst for item in sublist]
    trainr = [item for sublist in trlst for item in sublist]
    print('Done with making lists!')

    #for ind in range(len(flatr)):
    for ind in range(3):
        flatr[ind] = onehot(flatr[ind], biglist)
        featstr = '\t'.join(str(x) for x in flatr[ind][3])
        if frm == 'y':
            c, h, o, n, p, s, mim = fp.parse_form(flatr[ind][4])
            h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
            mad, amd, rmd = fp.mass_defects(mim)
            ffeatstr = '\t'.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
             h2o, c2p, c2n, c2s, amd, mad, rmd])
            testx.write(f'{featstr}\t{ffeatstr}\n')
        else: 
            testx.write(f'{featstr}\n')
        testy.write(f'{flatr[ind][1]}\n')
    testx.close()
    testy.close()
    print('Done with test!')

    #for ind in range(len(trainr)):
    for ind in range(3):
        trainr[ind] = onehot(flatr[ind], biglist)
        featstr = '\t'.join(str(x) for x in flatr[ind][3])
        if frm == 'y':
            c, h, o, n, p, s, mim = fp.parse_form(flatr[ind][4])
            h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
            mad, amd, rmd = fp.mass_defects(mim)
            ffeatstr = '\t'.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
             h2o, c2p, c2n, c2s, amd, mad, rmd])
            trainx.write(f'{featstr}\t{ffeatstr}\n')
        else:
            trainx.write(f'{featstr}\n')
        trainy.write(f'{trainr[ind][1]}\n')
    trainx.close()
    trainy.close()
    print('Done with train!')

def main():
    if len(sys.argv) == 1 or '-h' in sys.argv:
        show_help()
        sys.exit()
    
    try:
        infil, binw, outf, per, r, a, fr = argparsr(sys.argv)

    except: 
        print('Error reading arguments! Quitting!')
        show_help()
        sys.exit()

    mgffd = parse_mgf(open(infil, 'r'), fr)

    if a == 'y':
        make_arff(binw, mgffd, outf, per, fr)
    if r == 'y':
        make_irf(binw, mgffd, outf, per, fr)

    print('All done! Goodbye!')

if __name__ == '__main__':
    main()
