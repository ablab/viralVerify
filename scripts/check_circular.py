def check_circular(file, name):
   import fastaparser
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
            output.write(contig[1]+contig[1][kval:5000]+"\n")

   return(circular_contigs)


def main():
    import sys
    check_circular(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()

