import sys

print("INP1: Input file")
print("INP2: Index for class")
print("INP3: Name of positive class for binary classification")
file1=open(sys.argv[1],'r')
index=int(sys.argv[2])
cname=sys.argv[3]
line1=file1.readline()

#Make ARFF file
out1=open(sys.argv[1]+".%s.arff"%cname,'w')
out1x=open(sys.argv[1]+".%s.arff.instanceIDs"%cname,'w')
biglist=range(50,2500,1)
out1.write('@relation\tMSMS_binary_class\n\n')
for item in biglist:
    out1.write('@attribute\t%s\tnumeric\n'%item)
out1.write('@attribute\tclass\t{%s,not_%s}\n\n'%(cname,cname))
out1.write('@data\n')
pos=0; neg=0; m=0
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        
        em=int(float(tab1[7]))
        #print tab1, em
        #Get the class name
        cat=tab1[index]
        if cat==cname:
          ncat=cname
          pos+=1
        else:
          ncat=('not_%s'%cname)
          neg+=1

        #Get peaks
        peaks=tab1[8].split(',')
        plist=sorted([int(float(peak)) for peak in peaks])
        plist.append(em)        

        #One hot encoding
        encoded=[]
        for item in biglist:
            if item in plist:
                encoded.append('1')
            else:
                encoded.append('0')
        hot=','.join(encoded)
        out1.write('%s,%s\n'%(hot,ncat))
        out1x.write('%s\t%s\n'%(m,tab1[0]))


        m+=1
        if m%1000==0:
            print m, pos, neg
          
        
    line1=file1.readline()
file1.close(); out1.close(); out1x.close()

file1=open(sys.argv[1]+".%s.arff"%cname,'r')
out1=open(sys.argv[1]+".%s.arff.TEST.arff"%cname,'w')
out2=open(sys.argv[1]+".%s.arff.MAIN.arff"%cname,'w')
line1=file1.readline()
m=0; n=0
while line1:
    if line1.startswith('#') or line1.startswith('@') or line1=='':
        out1.write(line1)
        out2.write(line1)
    elif line1.strip().split()==[]:
        out1.write(line1)
        out2.write(line1)
    else:
        sp=line1.strip().split(',')
        if sp[-1]==cname:
            if m<20:                
                out1.write('%s,not_%s\n'%(','.join(sp[0:-1]),cname))
                m+=1
            else:
                out2.write(line1)
        else:
            if n<20:                
                out1.write('%s,%s\n'%(','.join(sp[0:-1]),cname))
                n+=1
            else:
                out2.write(line1)
    line1=file1.readline()
file1.close(); out1.close(); out2.close()
            

print "Positive: ", pos
print "Negative: ", neg
print "Total out: ", m
print "Done!"
          
          
        

