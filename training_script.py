# Training script for viralVerify
# Input - fasta files with viral, chromosomal and plasmid (optionally) training sequences

import sys
import argparse
import os, errno
from viralverify import check_circular

def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Training script for viralVerify")
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-v', required = True, help='File with viral nucleotide sequences in fasta format')
    required_args.add_argument('-nv', required = True, help='File with chromosomal nucleotide sequences in fasta format') 
    required_args.add_argument('--hmm', help='Path to HMM database') 
    required_args.add_argument('-o', required = True, help='Output directory')
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('-p', help='File with plasmid nucleotide sequences in fasta format')  
    optional_args.add_argument('-t', help='Number of threads')   
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()



args = parse_args(sys.argv[1:])
outdir = args.o
dirname = os.path.dirname(__file__)
input_files = {}

for i in ["v","nv","p"]:
  if (i in args) and getattr(args,i) != False:
    base_i = os.path.basename(getattr(args,i))
    input_files[i] = os.path.join(outdir, os.path.splitext(base_i)[0])

print (input_files)

try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise    

if args.hmm:
        hmm = args.hmm
else:
        print ("No HMM database provided") 
        exit(1) 
    
if args.t:
      threads = str(args.t)
else:
      threads = str(20)


# Run gene prediction

print ("Gene prediction...")

for i in input_files:
  res = os.system ("prodigal -p meta -c -i " + getattr(args,i) + " -a "+input_files[i]+"_proteins.fa -o "+input_files[i]+"_genes.fa 2>"+input_files[i]+"_prodigal.log" )
  if res != 0:
    print ("Prodigal run failed")
    exit(1)    

# HMM search

print ("HMM domains prediction...")

for i in input_files:
  res = os.system ("hmmsearch  --noali --cut_nc  -o "+input_files[i]+"_out_pfam --domtblout "+input_files[i]+"_domtblout --cpu "+ threads + " " + hmm + " "+input_files[i]+"_proteins.fa")
  if res != 0:
    print ("hmmsearch run failed")
    exit(1)


    

from operator import itemgetter


pfam={}
with open("/Nancy/mrayko/db/pfam/pfam_names.list", 'r') as infile:
    for line in infile:
      pfam[line.strip()]=[0,0,0]



#filter by e-value 
if "pl" in input_files:
  with open(input_files["pl"]+"_domtblout") as f:
    tblout_pl = f.read().splitlines()
    tblout_pl = [i.split() for i in tblout_pl]
    pl_count = 0
    for i in tblout_pl:
      if i[0][0]!= '#':
        if float(i[6]) <= 1e-06:
          pl_count+=1
          a = i[3]
          pfam[a][0] += 1

else:
  pl_count = 0


with open(input_files["nv"]+"_domtblout") as f1:
    tblout_chr = f1.read().splitlines()
    tblout_chr = [i.split() for i in tblout_chr]
    chr_count = 0
    for i in tblout_chr:
      if i[0][0]!= '#':
        if float(i[6]) <= 1e-06:
          chr_count+=1
          a = i[3]
          pfam[a][1] += 1



with open(input_files["v"]+"_domtblout") as f2:
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
table=[]
for key,value in pfam.items():
#  if value[0] >= 10 or value[1] >= 10:
    a = (value[0] + pseudo)/float(pl_count+len(pfam)*pseudo)
    b = (value[1] + pseudo)/float(chr_count+len(pfam)*pseudo)
    c = (value[2] + pseudo)/float(vir_count+len(pfam)*pseudo)
    d = (value[0] + value[1] + pseudo)/float(pl_count + chr_count +len(pfam)*pseudo) # for all non-viral

    table.append([key, value[0], value[1], value[2], a, b, c, d])

for i in sorted(table, key=itemgetter(5), reverse = True):
 # if i[1] >= 10 or i[2] >= 10 or i[3] >= 10:
    i = [str(j) for j in i]
    print ('\t'.join(i))
