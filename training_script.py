# Training script for PlasmidVerify
# Input - results of the hmmsearch (domtblout files) for plasmid, viral and chromosomal training dataset

import sys


plasmid_hmm = sys.argv[1]
virus_hmm = sys.argv[2]
chromosome_hmm = sys.argv[3]


from operator import itemgetter


pfam={}
with open("/Nancy/mrayko/db/pfam/pfam_names.list", 'r') as infile:
    for line in infile:
      pfam[line.strip()]=[0,0,0]



#filter by e-value 

with open(plasmid_hmm) as f:
    tblout_pl = f.read().splitlines()
    tblout_pl = [i.split() for i in tblout_pl]
    pl_count = 0
    for i in tblout_pl:
      if i[0][0]!= '#':
        if float(i[6]) <= 1e-06:
          pl_count+=1
          a = i[3]
          pfam[a][0] += 1

with open(chromsome_hmm) as f1:
    tblout_chr = f1.read().splitlines()
    tblout_chr = [i.split() for i in tblout_chr]
    chr_count = 0
    for i in tblout_chr:
      if i[0][0]!= '#':
        if float(i[6]) <= 1e-06:
          chr_count+=1
          a = i[3]
          pfam[a][1] += 1



with open(virus_hmm) as f2:
    tblout_vir = f2.read().splitlines()
    tblout_vir = [i.split() for i in tblout_vir]
    vir_count = 0
    for i in tblout_vir:
      if i[0][0]!= '#':
        if float(i[6]) <= 1e-06:
          vir_count+=1
          a = i[3]
          pfam[a][2] += 1


# add pseudocounts
pseudo=0.01
print (len(pfam)*pseudo)
table=[]
for key,value in pfam.items():
#  if value[0] >= 10 or value[1] >= 10:
    a = (value[0] + pseudo)/float(pl_count+len(pfam)*pseudo)
    b = (value[1] + pseudo)/float(chr_count+len(pfam)*pseudo)
    c = (value[2] + pseudo)/float(vir_count+len(pfam)*pseudo)
    d = (value[0] + value[1] + pseudo)/float(pl_count + chr_count +len(pfam)*pseudo) # for all non-viral

    table.append([key, value[0], value[1], value[2], a, b, c, d])

for i in sorted(table, key=itemgetter(5), reverse = True):
  if i[1] >= 10 or i[2] >= 10 or i[3] >= 10:
    i = [str(j) for j in i]
    print ('\t'.join(i))
