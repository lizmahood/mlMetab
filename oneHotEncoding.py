import sys, numpy as np, random
from collections import Counter

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

    return(inpf, binw, outp, perc)

def show_help():
    print('ARGS:\n'\
        '-inp full path to input file\n'\
        '-output full path to desired output file\n'\
        '-perct what percentage (as decimal) do you want for the test set?\n'\
        '-binw what is your preferred binwidth?')

def parse_mgf(inp):

    ##inp is a file object of the input mgf file
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
                mgfd[clas].append([name, clas, mz, peaks])
            else:
                mgfd[clas] = [[name, clas, mz, peaks]]

            ##resetting variables for next entry
            nam = inpl[6:].strip()
            clas = nam.split()[0].replace('[','').replace(']','')
            name = nam.replace(' ', '_').replace(':', '-').replace(';','').replace('}','').replace(',', '-')
        elif inpl.startswith('MW:'):
            mz = float(inpl.strip().split(': ')[1])
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

def make_arff(bnw, mgfd, out, p):
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
    for ind in range(len(flatr)):
        flatr[ind] = onehot(flatr[ind], biglist)
        featstr = ','.join(str(x) for x in flatr[ind][3])
        test.write(f'{flatr[ind][0]},{featstr},{flatr[ind][1]}\n')
    test.close()
    print('Done with test!')

    ##now training
    for ind in range(len(trainr)):
        trainr[ind] = onehot(trainr[ind], biglist)
        featstr = ','.join(str(x) for x in trainr[ind][3])
        train.write(f'{trainr[ind][0]},{featstr},{trainr[ind][1]}\n')
    train.close()
    print('Done with train!')
    
def main():
    if len(sys.argv) == 1 or '-h' in sys.argv:
        show_help()
        sys.exit()
    
    try:
        infil, binw, outf, per = argparsr(sys.argv)

    except: 
        print('Error reading arguments! Quitting!')
        show_help()
        sys.exit()

    mgffd = parse_mgf(open(infil, 'r'))
    make_arff(binw, mgffd, outf, per)
    print('All done! Goodbye!')

if __name__ == '__main__':
    main()
