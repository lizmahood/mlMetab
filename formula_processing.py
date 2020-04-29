from __future__ import division
import re, math

def returnval(list1):
     if len(list1)>0:
          return sum(list1)
     else:
          return 0.0

def getrat(a, b):
     '''
     Function to find ratios of elements. Avoids float division by 0. If the
     denominator is 0, this function will return the int 0.
     a and b can be floats or ints.
     Liz Mahood 12/12/17
     '''

     if b != 0 and b != 0.0:
          newval = a/b
     else:
          newval = 0

     return newval

def parse_form(frm):

    '''
    Function for finding numbers of elements in a formula
    and the rounded MIM of the formula.
    RETURNS: tuple with all CHNOPS values and mim
    '''
    good_els = {'C':12.000, 'H': 1.0078, 'O': 15.995, 
            'N': 14.003,  'P': 30.973, 'S': 31.972}

    elnums = re.findall('[A-Z][^A-Z]*', frm)

    numofC = 0; numofH = 0; numofN = 0
    numofO = 0; numofP = 0; numofS = 0

    for element in elnums:
        brokenel = re.findall('\d+|\D+', element)
        numofD = 0
        if brokenel[0] == 'D':
            try: 
                numofD = int(brokenel[1])
            except:
                print (frm)
            ##adds the # of Ds in this formula to # of Hs
            numofH = numofH + int(numofD)
        elif brokenel[0] == 'C':
            if len(brokenel) == 2:
                numofC = int(brokenel[1])
            elif len(brokenel) == 1:
                numofC = int(1)
        elif brokenel[0] == 'H':
            if len(brokenel) == 2:
                numofH = int(brokenel[1])
            elif len(brokenel) == 1:     
                numofH = int(1)
        elif brokenel[0] == 'O':
            if len(brokenel) == 2:
                numofO = int(brokenel[1])
            elif len(brokenel) == 1:
                numofO = int(1)
        elif brokenel[0] == 'N':
            if len(brokenel) == 2:
                numofN = int(brokenel[1])
            elif len(brokenel) == 1:
                numofN = int(1)
        elif brokenel[0] == 'P':
            if len(brokenel) == 2:
                numofP = int(brokenel[1])
            elif len(brokenel) == 1:
                numofP = int(1)
        elif brokenel[0] == 'S':
            if len(brokenel) == 2:
                numofS = int(brokenel[1])
            elif len(brokenel) == 1:
                numofS = int(1)
        else: continue
    
    mim = sum([numofC * good_els['C'], numofH * good_els['H'],
    numofO * good_els['O'], numofN * good_els['N'],
    numofP * good_els['P'], numofS * good_els['S']])
    mim = round(mim, 4)

    return(numofC, numofH, numofO, numofN, numofP, numofS, mim)

def mass_defects(mim):
    '''
    MIM is a float.
    RETURNS: tuple with mad, amd, rmd
    '''
    ##number after decimal
    mad = round(mim%1, 4)

    ##absolute mass defect
    roundmim = round(mim, 0)
    if round((roundmim- mim),4) >= 0:
        amd = round((roundmim - mim),4)
    elif round((roundmim - mim),4) < 0:
        amd = round(mad, 4)

    ##relative mass defect
    if mim != 0:
        rmd = int(round(((amd/mim) * 1e6),0))

    return(mad, amd, rmd)

def make_arff_header():
    ##RETURNS: multi-line string
    ##To be used as the arff header
    ##explaining the names and types of features

    feats = ['NewMIM', 'C', 'H', 'N', 'O', 'P', 'S',
    'CtoO', 'HtoC', 'HtoO', 'CtoP', 'CtoN', 'CtoS',
    'AMD', 'ADMD', 'RMD']
    arfstr = ''

    for feat in feats:
        arfstr += f'@attribute\t{feat}\tnumeric\n'

    return arfstr

def make_irf_header():
    feats = ['NewMIM', 'C', 'H', 'N', 'O', 'P', 'S',
    'CtoO', 'HtoC', 'HtoO', 'CtoP', 'CtoN', 'CtoS',
    'AMD', 'ADMD', 'RMD']

    ostr = '\t'.join(feats)
    return ostr

def get_allrat(c,h,o,n,p,s):
    '''
    Calcluates all 6 of our element ratios
    Returns them as floats
    '''
    h2c=getrat(h,c); h2o= getrat(h,o); c2o = getrat(c,o)
    c2n = getrat(c,n); c2s = getrat(c,s); c2p = getrat(c,p)

    return(h2c, h2o, c2o, c2n, c2s, c2p) 