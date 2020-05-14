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
        elif argl[ar] == '-chal':
            chal = argl[ar + 1]
        elif argl[ar] == '-classes':
            classes = argl[ar + 1]
    return(inpf, binw, outp, perc, makr, maka, frm, chal, classes)

def show_help():
    print('ARGS:\n'\
        '-inp full path to input file\n'\
        '-output full path to desired output file\n'\
        '-perct what percentage (as decimal) do you want for the test set?\n'\
        '-makeR do you want to make iRF input? y OR n\n'\
        '-makeA do you want to make WEKA input? y OR n\n'\
        '-classes FOR CHALLENGE, what classes were made in the training arff?\n'\
            'Leave blank if you don\'t want a challenge arff\n'\
        '-forms do you want to include formula features as well?\n'\
        '-chal do you want to make a challenge arff? n OR path to input\n'\
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
                    mgfd[clas] = [[name, clas, mz, peaks, form]]
                else:
                    mgfd[clas] = [[name, clas, mz, peaks]]

            ##resetting variables for next entry
            nam = inpl[6:].strip()
            clas = nam.split()[0].replace('[','').replace(']','')
            name = nam.replace(' ', '_').replace(':', '-').replace(';','').replace('}','').replace(',', '-')

        elif inpl.startswith('MW:'):
            mz = float(inpl.strip().split(': ')[1])
        elif inpl.startswith('Comment'):
            form = inpl.strip().split(';')[4].strip()
        elif inpl.startswith('Num Peaks:'):
            peaks = []
            inpl =inp.readline()
            while inpl[:1].isdigit():
                ##we are just logging mz, not abund
                peaks.append(float(inpl.strip().split()[0]))
                inpl = inp.readline()
    if frm == 'y':
        mgfd[clas].append([name, clas, mz, peaks, form])
    else: 
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

def make_outs(bnw, mgfd, out, p, frm, a, r):
    ##mgfd is the dict containing output of parse_mgf
    ##bnw is float or int
    ##out is output file name header, string
    ##p is a float, percent of file to reserve for test
    ## frm is y or n, strings
    ## a is y or n, strings
    ## r is y or n, strings
    ##OUTPUTS: arff and / or irf files with peaks as
    ##one-hot encoded features

    biglist=np.arange(50,2500,float(bnw))

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

    if a == 'y':
        #Initialize ARFF file
        train=open(out + 'training.arff','w')
        test=open(out + 'test.arff','w')

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

    if r == 'y':
        trainx = open(out + 'xtrain.tab', 'w')
        testx = open(out + 'xtest.tab', 'w')
        trainy = open(out + 'ytrain.tab', 'w')
        testy = open(out + 'ytest.tab', 'w')

        ##writing out feature names for x-files
        news = ''
        for i in range(len(biglist)-1):
            news += f'{biglist[i]}-{biglist[i+1]}\t'
    
        if frm == 'y':
            news += fp.make_irf_header()

        trainx.write(f'{news}\n')
        testx.write(f'{news}\n')

    ##doing one-hot now and writing out
    ##first for test arff instances
    for ind in range(len(flatr)):
        flatr[ind] = onehot(flatr[ind], biglist)
        if a == 'y':
            featstra = ','.join(str(x) for x in flatr[ind][3])
        if r == 'y':
            featstrr = '\t'.join(str(x) for x in flatr[ind][3])
        if frm == 'y':
            #print(flatr[ind])
            try:
                c, h, o, n, p, s, mim = fp.parse_form(flatr[ind][4])
            except IndexError:
                print(flatr[ind])
            h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
            mad, amd, rmd = fp.mass_defects(mim)
            if a == 'y':
                ffeatstra = ','.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
                h2o, c2p, c2n, c2s, amd, mad, rmd])
                test.write(f'{flatr[ind][0]},{featstra},{ffeatstra},{flatr[ind][1]}\n')
            if r == 'y':
                ffeatstrr = '\t'.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
                h2o, c2p, c2n, c2s, amd, mad, rmd])
                testx.write(f'{featstrr}\t{ffeatstrr}\n')
        else:
            if a == 'y':
                test.write(f'{flatr[ind][0]},{featstra},{flatr[ind][1]}\n')
            if r == 'y':
                testx.write(f'{featstrr}\n')
        if r == 'y':
            testy.write(f'{flatr[ind][1]}\n')
    if a == 'y':
        test.close()
    if r == 'y':
        testx.close()
        testy.close()
    print('Done with test(s)!')

    ##now training
    for ind in range(len(trainr)):
        trainr[ind] = onehot(trainr[ind], biglist)
        if a == 'y':
            featstra = ','.join(str(x) for x in trainr[ind][3])
        if r == 'y':
            featstrr = '\t'.join(str(x) for x in trainr[ind][3]) 
        if frm == 'y':
            try:
                c, h, o, n, p, s, mim = fp.parse_form(trainr[ind][4])
            except IndexError:
                print(trainr[ind])
            h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
            mad, amd, rmd = fp.mass_defects(mim)
            if a == 'y':
                ffeatstra = ','.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
                h2o, c2p, c2n, c2s, amd, mad, rmd])
                train.write(f'{trainr[ind][0]},{featstra},{ffeatstra},{trainr[ind][1]}\n')
            if r == 'y':
                ffeatstrr = '\t'.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
                h2o, c2p, c2n, c2s, amd, mad, rmd])
                trainx.write(f'{featstrr}\t{ffeatstrr}\n')
        else: 
            if a == 'y':
                train.write(f'{trainr[ind][0]},{featstra},{trainr[ind][1]}\n')
            if r == 'y':
                trainx.write(f'{featstrr}\n')
        if r == 'y':
            trainy.write(f'{trainr[ind][1]}\n')
    if a == 'y':
        train.close()
    if r == 'y':
        trainx.close()
        trainy.close()
    print('Done with train(s)!')

def make_challenge_arff(inp, frm, out, binw, clss):
    ##inp is the path to the mgf file
    ##Frm is a string (y OR n)
    ##out is path to output name
    ##binw is int or float
    ##clss is large string of all classes used in training arff
    inp = open(inp, 'r', encoding = 'utf-8')
    chal = open(out + '_challenge.arff', 'w')
    biglist=np.arange(50,2500,float(binw))

    #Initialize ARFF file
    chal.write('@relation\tMSMS_multiclass\n\n')
    chal.write('@attribute\tID\tstring\n')

    
    for i in range(len(biglist)-1):
        chal.write(f'@attribute\t{biglist[i]}-{biglist[i+1]}\tnumeric\n')

    if frm == 'y':
        frmft = fp.make_arff_header()
        chal.write(f'{frmft}')

    ##adding line for all possible classes -- in training arff
    chal.write(f'@attribute\tclass\t{{{clss}}}\n\n')
    chal.write('@data\n')

    ##looping through challenge msp
    strt = True
    inpl = inp.readline()
    while inpl:
        olst = []
        if (inpl.startswith('Name') and strt == True) or (inpl.startswith('SCAN') and strt == True):
            nam = inpl[6:].strip()
            name = nam.replace(' ', '_').replace(':', '-').replace(';', '').replace('}','')\
                .replace(',','-').replace('γ', 'Gam').replace('⁴-⁷','').replace('"', '')\
                    .replace('⁴-⁹','').replace('(−)','').replace('²-⁷','')\
                        .replace('³-⁵','').replace('′','').replace('\'','').replace('{','')
            clas = '?'
            strt = False
            if inpl.startswith('SCAN'):
                form = inpl.strip().split('|')[2]

        elif (inpl.startswith('Name') and strt != True) or (inpl.startswith('SCAN') and strt != True):
            if frm == 'y':
                olst = [name, clas, mz, peaks, form]
            else:
                olst = [name, clas, mz, peaks]

            nlst = onehot(olst, biglist)
            featstra = ','.join(str(x) for x in nlst[3])
            if frm == 'y':
                c, h, o, n, p, s, mim = fp.parse_form(nlst[4])
                h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
                mad, amd, rmd = fp.mass_defects(mim)
                ffeatstra = ','.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
                h2o, c2p, c2n, c2s, amd, mad, rmd])
                try:
                    chal.write(f'{nn},{featstra},{ffeatstra},{nlst[1]}\n')  
                except:
                    encc = nlst[0].encode('ASCII', 'ignore')
                    encd = encc.decode()
                    chal.write(f'{encd},{featstra},{ffeatstra},{nlst[1]}\n')
          
            else:
                chal.write(f'{nlst[0]},{featstra},{nlst[1]}\n')

            ##resetting variables for next entry
            nam = inpl[6:].strip()
            clas = '?'
            name = nam.replace(' ', '_').replace(':', '-').replace(';', '').replace('}','')\
                .replace(',','-').replace('γ', 'Gam').replace('⁴-⁷','').replace('"', '')\
                    .replace('⁴-⁹','').replace('(−)','').replace('²-⁷','')\
                        .replace('³-⁵','').replace('′','').replace('\'','').replace('{','')
            if inpl.startswith('SCAN'):
                form = inpl.strip().split('|')[2]
        elif inpl.startswith('PrecursorMZ:'):
            mz = float(inpl.strip().split(': ')[1])
        elif inpl.startswith('PEPMASS'):
            mz = float(inpl.strip().split('=')[1])
        elif inpl.startswith('Formula'):
            form = inpl.strip().split(': ')[1]
        elif inpl.startswith('Num Peaks:'):
            peaks = []
            inpl =inp.readline()
            while inpl[:1].isdigit():
                ##we are just logging mz, not abund
                peaks.append(float(inpl.strip().split()[0]))
                inpl = inp.readline()
        elif inpl.startswith('ION='):
            peaks = []
            inpl = inp.readline()
            while inpl[1].isdigit():
                peaks.append(float(inpl.strip().split()[0]))
                inpl = inp.readline()
        inpl = inp.readline()
    if frm == 'y':
        olst = [name, clas, mz, peaks, form]
    else:
        olst = [name, clas, mz, peaks]

    nlst = onehot(olst, biglist)
    featstra = ','.join(str(x) for x in nlst[3])
    if frm == 'y':
        c, h, o, n, p, s, mim = fp.parse_form(nlst[4])
        h2c, h2o, c2o, c2n, c2s, c2p = fp.get_allrat(c,h,o,n,p,s)
        mad, amd, rmd = fp.mass_defects(mim)
        ffeatstra = ','.join(str(x) for x in [mim, c, h, n, o, p, s, c2o, h2c,
        h2o, c2p, c2n, c2s, amd, mad, rmd])
        chal.write(f'{nlst[0]},{featstra},{ffeatstra},{nlst[1]}\n')            
    else:
        chal.write(f'{nlst[0]},{featstra},{nlst[1]}\n')

def main():
    if len(sys.argv) == 1 or '-h' in sys.argv:
        show_help()
        sys.exit()
    
    try:
        infil, binw, outf, per, r, a, fr, chal, classes = argparsr(sys.argv)

    except: 
        print('Error reading arguments! Quitting!')
        show_help()
        sys.exit()

    if r != 'n' and a != 'n':
        mgffd = parse_mgf(open(infil, 'r'), fr)
        make_outs(binw, mgffd, outf, per, fr, a, r)
    
    if chal != 'n':
        make_challenge_arff(chal, fr, outf, binw, classes)

    print('All done! Goodbye!')

if __name__ == '__main__':
    main()
