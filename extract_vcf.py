# get vcf file for a specific region.
# inputs: vcffile.vcf #chr start end

import sys
import re
import string

if len(sys.argv) != 6:
    print "Usage: python extract_vcf.py [inputvcf] [#chr] [start] [end] [output]"
    sys.exit (1)
    
input = open(sys.argv[1], 'r')
chr = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])
output = open(sys.argv[5], 'w')

for l in input:
    l = l.strip()
    if l.startswith('#'):
        output.write(l+'\n')
        continue
    fields = l.split('\t')
    if fields[0] == 'chr'+str(chr): # depend on what you have in your vcf file, if there is no 'chr', then delete the 'chr'
        if int(fields[1]) >= start:
            if int(fields[1]) <= end:
                output.write(l+'\n')
                
output.close()
