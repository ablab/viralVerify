#!/usr/bin/env python
import os, errno
import sys
import argparse
import collections
from math import log
from math import exp
import csv
import operator
import time
import datetime
import fastaparser


def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="HMM-based plasmid verification script")
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-f', required = True, help='Input fasta file')
    required_args.add_argument('-o', required = True, help='Output directory')
    required_args.add_argument('--hmm', help='Path to Pfam-A HMM database')    
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--db', help='Run BLAST on input contigs with provided database')
    optional_args.add_argument('-t', help='Number of threads')   
    optional_args.add_argument('-thr', help='Detection threshold for classifier (default = 7)')   
    optional_args.add_argument('-p', action='store_true', help='Output predicted plasmids separately')   
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()



def check_circular(file, name):
   contigs = fastaparser.read_fasta(file)
   count = []
   circular_contigs = {}
   input_fasta = name + "_input_with_circ.fasta"

   with open(input_fasta, 'w') as output:

    for contig in contigs:
       circular_contigs[contig[0].split(" ")[0][1:]] = [len(contig[1]), "-"]
       for kval in range (200,50, -1):
           if kval >= len(contig[1]) or len(contig[1]) < 500:
               continue
           start = contig[1][:kval]
           end = contig[1][-kval:]

           if start == end: 
               circular_contigs[contig[0].split(" ")[0][1:]] = [len(contig[1]), "+"]
               break
       

       output.write(contig[0]+"\n")
       if circular_contigs[contig[0].split(" ")[0][1:]][1] == "-":
          output.write(contig[1]+"\n")
       elif circular_contigs[contig[0].split(" ")[0][1:]][1] == "+":
        #  if contig[0].split(" ")[0][1:] in circ_set:
          #  output.write(contig[1]+contig[1][:5000]+"\n")
         # else:
            output.write(contig[1]+contig[1][kval:5000]+"\n")

   return(circular_contigs)
    


def get_table_from_tblout(tblout_pfam):
    with open(tblout_pfam, "r") as infile:
        tblout_pfam=infile.readlines()
   
    tblout_pfam = [i.split() for i in tblout_pfam[3:-10]]
    for i in tblout_pfam:
        i[13] = float(i[13])

    tblout_pfam.sort(key = operator.itemgetter(0, 13,17), reverse = True)

    top_genes={}
    for i in tblout_pfam:
        if i[0] not in top_genes:
            top_genes[i[0]] = [[i[3],float(i[13]),float(i[17]),float(i[18])]]
        else:
            for j in top_genes[i[0]]:
                start_i, end_i, start_j, end_j = float(i[17]), float(i[18]), float(j[2]), float(j[3])
                 
                if not ((end_i <= start_j) or (start_i >= end_j)):
                    break
                else: 
                    top_genes[i[0]].append([i[3],float(i[13]),start_i,end_i])
                    break


    contigs = collections.OrderedDict()
    for i in top_genes:
        name = i.rsplit("_", 1)[0]
        if name not in contigs:
            contigs[name] = set()
            for i in top_genes[i]:
                contigs[name].add(i[0])
        else:
            for i in top_genes[i]:
                contigs[name].add(i[0])

    out = []
    for key, value in contigs.items():
        out+=[str(key) + " "  +  " ".join(value)]

    return out


def naive_bayes(input_list):
    unc_score = threshold
    tr=os.path.dirname(os.path.abspath(__file__)) + "/classifier_table.txt"

    with open(tr, 'r') as infile:
        table=infile.readlines()
        table = [i.split() for i in table]

# hmm dictionary - for each HMM store plasmid, chromosomal, viral and pl+chr frequency 
    hmm_dict = {}
    for i in table:
        hmm_dict[i[0]] = [float(i[4]),float(i[5]),float(i[6]), float(i[7])]

# Calculate probabilities for each element of input list
    out_list=[]
    for i in input_list:
        chrom, plasm,vir, comb, chrom_log, plasm_log, vir_log, comb_log = 1, 1, 1, 1, 0, 0, 0, 0
        for j in i.split():
          if j in hmm_dict.keys(): 
                plasm=plasm*hmm_dict[j][0]
                plasm_log=plasm_log+log(hmm_dict[j][0])
                chrom=chrom*hmm_dict[j][1]
                chrom_log=chrom_log+log(hmm_dict[j][1])
                vir=vir*hmm_dict[j][2]
                vir_log=vir_log+log(hmm_dict[j][2])
                comb=comb*hmm_dict[j][3]
                comb_log=comb_log+log(hmm_dict[j][3])
        if (vir_log - comb_log) > unc_score: 
             out_list.append(["Virus", vir_log, comb_log, "{0:.2f}".format(vir_log - comb_log)])
        elif (vir_log - comb_log) > (-1)*unc_score:
          if len(i.split()) > 2:
             out_list.append(["Uncertain - viral or bacterial", vir_log, comb_log, "{0:.2f}".format(vir_log - comb_log)])
          else:
             out_list.append(["Uncertain - too short", vir_log, comb_log, "{0:.2f}".format(vir_log - comb_log)])

        elif (plasm_log - chrom_log) > unc_score:
             out_list.append(["Plasmid", plasm_log, chrom_log, "{0:.2f}".format(plasm_log - chrom_log)])
        elif (chrom_log - plasm_log) > unc_score:
             out_list.append(["Chromosome", vir_log, comb_log, "{0:.2f}".format(plasm_log - chrom_log)])
        else: 
             out_list.append(["Uncertain - plasmid or chromosomal", vir_log, comb_log, "{0:.2f}".format(plasm_log - chrom_log)])



        
    return out_list 



def main():

    args = parse_args(sys.argv[1:])
    base = os.path.basename(args.f)
    name_file = os.path.splitext(base)[0]
    dirname = os.path.dirname(__file__)
    outdir = args.o
    
    try:
        os.makedirs(outdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    
    name = os.path.join(outdir, name_file)
    
    ids = []
    with open(args.f, "r") as ins:
        for line in ins:
            if line[0]==">":
                ids.append(line.split()[0][1:])
    
    if args.hmm:
        hmm = args.hmm
    else:
        print ("No HMM database provided") 
        exit(1)    
    
    
    if args.db:
        from parse_blast_xml import parser
        blastdb = args.db
    
    if args.t:
        threads = str(args.t)
    else:
        threads = str(20)

    if args.thr:
        threshold = int(args.thr)
    else:
        threshold = 7 
    

    # Check for circular:   
    contig_len_circ = check_circular(args.f, name)
    infile_circ = name + "_input_with_circ.fasta"

    # Run gene prediction
    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

    print ("Gene prediction...")
    res = os.system ("prodigal -p meta -c -i " + infile_circ + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
    if res != 0:
        print ("Prodigal run failed")
        exit(1)    

    # Filter genes predicted over the end of the contig
  
    proteins = fastaparser.read_fasta(name+"_proteins.fa")
    with open(name+"_proteins_circ.fa", 'w') as protein_output:
      for i in proteins:
        contig_name = i[0].split()[0].rsplit("_",1)[0][1:]
        gene_start = i[0].split("#")[1]
        if int(gene_start.strip()) < int((contig_len_circ[contig_name][0])):
          protein_output.write(i[0]+"\n")
          protein_output.write(i[1]+"\n")

   # HMM search

    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
    print ("HMM domains prediction...")
    res = os.system ("hmmsearch  --noali --cut_nc  -o "+name+"_out_pfam --domtblout "+name+"_domtblout --cpu "+ threads + " " + hmm + " "+name+"_proteins_circ.fa")
    if res != 0:
        print ("hmmsearch run failed")
        exit(1)  

    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
   
    print ("Parsing...")
    tblout_pfam= name + "_domtblout" 


    feature_table = get_table_from_tblout(tblout_pfam) 
    feature_table = [i.strip().split(' ', 1) for i in feature_table]
    
    with open(name + '_feature_table.txt', 'w') as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows(feature_table)
    
    
    feature_table_names=[]
    feature_table_genes=[]
    for i in feature_table:
          feature_table_names.append(i[0])
          feature_table_genes.append(i[1])
    
    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))  
    
    print ("Classification...")
    t=feature_table_genes
    k = naive_bayes(t)

    names_result={}
    for i in range (0,len(k)):
      names_result[feature_table_names[i]] = [k[i][0],k[i][3],feature_table_genes[i]]

    


    if args.db:
        #run blast
        print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
        print ("Running BLAST...")
    
        os.system ("blastn  -query " + args.f + " -db " + blastdb + " -evalue 0.0001 -outfmt 5 -out "+name+".xml -num_threads "+threads+" -num_alignments 50")
        print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
        print ("Parsing BLAST")
        parser(name+".xml", outdir)
        print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

    
        #### add blast results
        plasmids= [line.strip().split("\t") for line in open(name + "_plasmid.names")]
    
        plasmids_list={}
        for i in range(0, len(plasmids)-1):
          if len(plasmids[i])==1:
            plasmids_list[plasmids[i][0].split()[0]] = [float(plasmids[i+1][1].split(":")[1]), float(plasmids[i+1][2].split(":")[1]), plasmids[i+1][0]]
    
        chrom= [line.rstrip().split("\t") for line in open(name + "_chromosome.names")]
        chrom_list={}
        for i in range(0, len(chrom)-1):
          if len(chrom[i])==1:
            chrom_list[chrom[i][0].split()[0]] = [float(chrom[i+1][1].split(":")[1]), float(chrom[i+1][2].split(":")[1]), chrom[i+1][0]]
    
    
        vir= [line.rstrip().split("\t") for line in open(name + "_viruses.names")]
        vir_list={}
        for i in range(0, len(vir)-1):
          if len(vir[i])==1:
                vir_list[vir[i][0].split()[0]] = [float(vir[i+1][1].split(":")[1]), float(vir[i+1][2].split(":")[1]), vir[i+1][0]]
    
        nos= [line.rstrip() for line in open(name + "_no_significant.names")]
        nos_list=[]
        for i in nos:
          if len(i.split())>0:
            nos_list.append(i.split()[0])
        nos_list = [i.strip().split()[0] for i in nos_list]
    
    
        other= [line.rstrip() for line in open(name + "_other.names")]
        other_list=[]
        for i in other:
          if len(i.split())>0:
            other_list.append(i.split()[0])
        other_list = [i.strip().split()[0] for i in other_list] 

    
    final_table=collections.OrderedDict()
    if args.db:
     for i in ids:
       if i in names_result:
            if names_result[i][0] == "Uncertain - too short":
              if (contig_len_circ[i][0] > 3000) or (contig_len_circ[i][1] == "+"):
                names_result[i][0] = "Uncertain - viral or bacterial"

            if i in plasmids_list:
              final_table[i] = [names_result[i][0], contig_len_circ[i][0], contig_len_circ[i][1], names_result[i][1], names_result[i][2], "Plasmid", round(plasmids_list[i][0],2), round(plasmids_list[i][1],2),plasmids_list[i][2]]
            if i in chrom_list:
              final_table[i] = [names_result[i][0], contig_len_circ[i][0], contig_len_circ[i][1], names_result[i][1], names_result[i][2], "Chromosome", chrom_list[i][0], chrom_list[i][1],chrom_list[i][2]]
            if i in vir_list:
              final_table[i] = [names_result[i][0], contig_len_circ[i][0], contig_len_circ[i][1], names_result[i][1], names_result[i][2], "Virus", vir_list[i][0], vir_list[i][1],vir_list[i][2]]
            if i in nos_list:
              final_table[i] = [names_result[i][0], contig_len_circ[i][0], contig_len_circ[i][1], names_result[i][1], names_result[i][2], "Non-significant"]
            if i in other_list:
              final_table[i] = [names_result[i][0], contig_len_circ[i][0], contig_len_circ[i][1], names_result[i][1], names_result[i][2], "Other", ]
    
       else: 
            if (contig_len_circ[i][0] > 3000) or (contig_len_circ[i][1] == "+"):
              names_result[i] = "Uncertain - viral or bacterial"
            else:
              names_result[i] = "Uncertain - too short"

              
            if i in plasmids_list:
              final_table[i] = [names_result[i], contig_len_circ[i][0], contig_len_circ[i][1], "-", "-", "Plasmid", plasmids_list[i][0], plasmids_list[i][1],plasmids_list[i][2]]
            if i in chrom_list:
              final_table[i] = [names_result[i], contig_len_circ[i][0], contig_len_circ[i][1], "-", "-", "Chromosome", chrom_list[i][0], chrom_list[i][1],chrom_list[i][2]]
            if i in vir_list:
              final_table[i] = [names_result[i], contig_len_circ[i][0], contig_len_circ[i][1],"-", "-", "Virus", vir_list[i][0], vir_list[i][1],vir_list[i][2]]
            if i in nos_list:
              final_table[i] = [names_result[i], contig_len_circ[i][0], contig_len_circ[i][1], "-", "-", "Non-significant"]
            if i in other_list:
              final_table[i] = [names_result[i], contig_len_circ[i][0], contig_len_circ[i][1], "-", "-", "Other", other_list[i][0], other_list[i][1],other_list[i][2]]



    else:
     for i in ids: 
      if i in names_result:
        if names_result[i][0] == "Uncertain - too short":
          if (contig_len_circ[i][0] > 3000) or (contig_len_circ[i][1] == "+"):
            names_result[i][0] = "Uncertain - viral or bacterial"
        final_table[i] = [names_result[i][0], contig_len_circ[i][0], contig_len_circ[i][1], names_result[i][1],names_result[i][2]]
      else:
        if (contig_len_circ[i][0] > 3000) or (contig_len_circ[i][1] == "+"):
          final_table[i] = ["Uncertain - viral or bacterial", contig_len_circ[i][0], contig_len_circ[i][1], "-"]
        else:
          final_table[i] = ["Uncertain - too short", contig_len_circ[i][0], contig_len_circ[i][1], "-"]
    
    result_file = name + "_result_table.csv"
    with open(result_file, 'w') as output:
        writer = csv.writer(output, lineterminator='\n')
        for i in final_table:
          writer.writerow([i] + final_table[i])
    
    if not os.path.exists(outdir + "/Prediction_results_fasta/"):
        os.mkdir(outdir + "/Prediction_results_fasta/")


    with open (outdir + "/Prediction_results_fasta/" +  name_file + "_virus.fasta", "w") as vir_file:
      with open (outdir + "/Prediction_results_fasta/" +  name_file + "_plasmid.fasta", "w") as plasmid_file:
        with open (outdir + "/Prediction_results_fasta/" +  name_file + "_chromosome.fasta", "w") as chrom_file:
          with open (outdir + "/Prediction_results_fasta/" +  name_file + "_virus_uncertain.fasta", "w") as vc_file:
            with open (outdir + "/Prediction_results_fasta/" +  name_file + "_plasmid_uncertain.fasta", "w") as pc_file:
                contigs = fastaparser.read_fasta(args.f)
                for i in contigs:
                  contig_name = i[0].split(" ")[0][1:]
                  if final_table[contig_name][0] == "Virus":
                    vir_file.write(i[0]+"\n")
                    vir_file.write(i[1]+"\n")
                  elif final_table[contig_name][0] == "Chromosome":
                    chrom_file.write(i[0]+"\n")
                    chrom_file.write(i[1]+"\n")                    
                  elif final_table[contig_name][0] == "Plasmid":
                    if args.p:
                      plasmid_file.write(i[0]+"\n")
                      plasmid_file.write(i[1]+"\n")
                    else:
                      chrom_file.write(i[0]+"\n")
                      chrom_file.write(i[1]+"\n")    

                  elif final_table[contig_name][0] == "Uncertain - viral or bacterial":
                    vc_file.write(i[0]+"\n")
                    vc_file.write(i[1]+"\n")   

                  elif final_table[contig_name][0] == "Uncertain - plasmid or chromosomal":
                    if args.p:
                      pc_file.write(i[0]+"\n")
                      pc_file.write(i[1]+"\n")
                    else:
                      chrom_file.write(i[0]+"\n")
                      chrom_file.write(i[1]+"\n")  

                if not args.p:
                    os.remove(outdir + "/Prediction_results_fasta/" +  name_file + "_plasmid.fasta")  
                    os.remove(outdir + "/Prediction_results_fasta/" +  name_file + "_plasmid_uncertain.fasta")  



    print ("Done!")
    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))  
    print ("Verification results can be found in " + os.path.abspath(result_file))
    
    

if __name__ == "__main__":
    main()

