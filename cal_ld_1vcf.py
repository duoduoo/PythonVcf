# calculate linkage (R2) from 1 vcf files, for two sets of variants selected by users.
# The second set linked to the first set.
# The vcf files should be from 1000 genome phase 3.
# only consider 0,1 match, output alternative allele variations, and for further consideration.
# exclude X chromosome!!!
# error is the output for variations with "2" alternative allele from both files.
# in phase 3, there are 2504 individuals.
#
##########################################################################################
#
# 
# SNP1  SNP2    R2
# rs1   rs2     ...
# rs1   rs3     ...
# ...
#
##########################################################################################

import sys
import re
import string

if len(sys.argv) != 6:
    print "Usage: python cal_ld_1vcf.py [SNP_vcf] [SNP1_ids] [SNP2_ids] [output] [error]"    #launch the program the the sequence file
    sys.exit (1)
    
vcf = open(sys.argv[1], 'r')
snp1 = open(sys.argv[2], 'r') # snp set 1 
snp2 = open(sys.argv[3], 'r') # snp set 2
output = open(sys.argv[4], 'w')
error = open(sys.argv[5], 'w')

output.write('SNP1\tSNP2\tR2\n')

snpSet1 = []
snpSet2 = []

# function counting 1 allele for a list of genotype
def AlleleCount(genotype_list):
    ac_dic = {}
    for genotype in genotype_list:
        if genotype not in ac_dic:
            ac_dic[genotype] = 1
        else:
            ac_dic[genotype] = ac_dic[genotype] + 1

    allele1_count = 0
    for element in ac_dic:
        if element == '1|1':
            allele1_count = allele1_count + ac_dic[element]*2
        elif element == '1|0':
            allele1_count = allele1_count + ac_dic[element]
        elif element == '0|1':
            allele1_count = allele1_count + ac_dic[element]

    if allele1_count <= 2504:
        return [allele1_count, 1] # minor allele is 1
    elif allele1_count > 2504:
        minor_allele = 5008-allele1_count
        return [minor_allele, 0] # minor allele is 0



# pos1 and pos2 here mean snp id.
for pos1 in snp1:
    pos1 = pos1.strip()
    snpSet1.append(pos1)

for pos2 in snp2:
    pos2 = pos2.strip()
    snpSet2.append(pos2)

snp_vcf = {}

for line in vcf:
    line = line.strip()
    if line.startswith('##'):
        continue
    if line.startswith('X'):
        continue
    if line.startswith('#'):
        id = []
        id = line.split('\t')
        continue
    gt = []
    gt = line.split('\t')
    if gt[2] in snpSet1:
        snp_vcf[gt[2]] = gt[9:2513]
    elif gt[2] in snpSet2:
        snp_vcf[gt[2]] = gt[9:2513]
    # this version, use SNP id, instead of position

for snp in snpSet1:
    snp_gt = snp_vcf[snp]

    snp_ac1 = AlleleCount(snp_gt)[0]
    snp_minor_allele = AlleleCount(snp_gt)[1]
    
    potential_ld_pos = []
    potential_ld_rev = []
    
    for sv in snpSet2:
        # ld matching 0|1 with 0|1 (positive)

        sv_ac1 = AlleleCount(snp_vcf[sv])[0]
        sv_minor_allele = AlleleCount(snp_vcf[sv])[1]

        if snp_minor_allele == 1 and sv_minor_allele == 1:
            match_1_count = 0
            for n in range(0, 2504):
                # match '1' in '1|0'
                if snp_gt[n][:1] == '1' and snp_vcf[sv][n][:1] == '1':
                    match_1_count = match_1_count + 1
                # match '0' in '1|0'
                if snp_gt[n][-1:] == '1' and snp_vcf[sv][n][-1:] == '1':
                    match_1_count = match_1_count + 1

        elif snp_minor_allele == 1 and sv_minor_allele == 0:
            match_1_count = 0
            for n in range(0, 2504):
                # match '1' in '1|0'
                if snp_gt[n][:1] == '1' and snp_vcf[sv][n][:1] == '0':
                    match_1_count = match_1_count + 1
                # match '0' in '1|0'
                if snp_gt[n][-1:] == '1' and snp_vcf[sv][n][-1:] == '0':
                    match_1_count = match_1_count + 1

        elif snp_minor_allele == 0 and sv_minor_allele == 1:
            match_1_count = 0
            for n in range(0, 2504):
                # match '1' in '1|0'
                if snp_gt[n][:1] == '0' and snp_vcf[sv][n][:1] == '1':
                    match_1_count = match_1_count + 1
                # match '0' in '1|0'
                if snp_gt[n][-1:] == '0' and snp_vcf[sv][n][-1:] == '1':
                    match_1_count = match_1_count + 1

        elif snp_minor_allele == 0 and sv_minor_allele == 0:
            match_1_count = 0
            for n in range(0, 2504):
                # match '1' in '1|0'
                if snp_gt[n][:1] == '0' and snp_vcf[sv][n][:1] == '0':
                    match_1_count = match_1_count + 1
                # match '0' in '1|0'
                if snp_gt[n][-1:] == '0' and snp_vcf[sv][n][-1:] == '0':
                    match_1_count = match_1_count + 1

        else:
            error.write('Something wrong!')

        fenmu = snp_ac1*(5008-snp_ac1)*sv_ac1*(5008-sv_ac1)
        fenzi = (match_1_count*5008-snp_ac1*sv_ac1)*(match_1_count*5008-snp_ac1*sv_ac1)

        r2 = float(fenzi) / float(fenmu)
        output.write(snp+'\t'+sv+'\t'+str(r2)+'\n')

    output.write('\n')


print "Erica is a genius!"    

output.close()
