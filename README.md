# viralVerify: viral contig verification tool

viralVerify classifies contigs (output of metaviralSPAdes or other assemblers) as viral, non-viral or uncertain, 
based on gene content. Also for non-viral contigs it can optionally provide plasmid/non-plasmid classification.

viralVerify predicts genes in the contig using Prodigal in the metagenomic mode, runs hmmsearch on the predicted proteins 
and classifies the contig as vrial or non-viral by applying the Naive Bayes classifier (NBC). 
For the set of predicted HMMs, viralVerify uses trained NBC to classify this set to be viral or chromosomal. 

### Requirements

viralVerify is a Python script, thus, installation is not required. However, it has the following dependencies:

* Python 2.7+,
* Prodigal (https://github.com/hyattpd/Prodigal, available via conda),
* hmmsearch (from the hmmer package, http://hmmer.org/download.html),
* [provided database of virus/chromosome-specific HMMs](http://data.cab.spbu.ru/index.php/s/RogrHJtxBRrY7ec/download?path=%2F&files=nbc_hmms.hmm.gz)
 or 
* Pfam-A database (you can download recent release here: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/).

To work properly, viralVerify require Prodigal and hmmsearch in your PATH environment variable.


### Optional BLAST verification

You can verify your output by BLAST to check if you found novel viruses or plasmids. In this case, you need to have blastn in your $PATH, Biopython installed, and provide a path to the local copy of nt database. 

## Usage 

    ./viralverify.py 
            -f Input fasta file
            -o output_directory 
            --hmm HMM  Path to HMM database

            Optional arguments:
            -h, --help  Show the help message and exit
            --db DB     Run BLAST on input contigs against provided db
            -t          Number of threads
            -p          Output predicted plasmidic contigs separately


Output file: comma-separated table <input_file>_result_table.csv

Output format: contig name, prediction result, log-likelihood ratio, list of predicted HMMs
  
Fasta files with prediction results can be found in Prediction_results_fasta folder
  
