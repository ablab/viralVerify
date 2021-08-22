# viralVerify: viral contig verification tool

**Version: 1.1**

viralVerify classifies contigs (output of metaviralSPAdes or other assemblers) as viral, non-viral or uncertain, 
based on gene content. Also for non-viral contigs it can optionally provide plasmid/non-plasmid classification.

viralVerify predicts genes in the contig using Prodigal in the metagenomic mode, runs hmmsearch on the predicted proteins 
and classifies the contig as vrial or non-viral by applying the Naive Bayes classifier (NBC). 
For the set of predicted HMMs, viralVerify uses trained NBC to classify this set to be viral or chromosomal. 

To improve results in the case of metagenomes with possible host contamination, we recommend users to filter out reads that align to the host genome prior to assembly.
Since viralVerify is based on gene classification, it can be used on contigs on any length, and short viruses can be detected as long as they contain a recognizable virus-specific gene. To help analyze the rapidly growing amount of novel data, we have added a script that allows users to construct their own training database from a set of viral, chromosomal and plasmid contigs, as well as custom HMM database

### Requirements

viralVerify is a Python script, thus, installation is not required. However, it has the following dependencies:

* Python 3.6+,
* Prodigal (https://github.com/hyattpd/Prodigal, available via conda),
* hmmsearch (from the hmmer package, http://hmmer.org/download.html),
* provided *decompressed* database of virus/chromosome-specific HMMs (https://figshare.com/s/f897d463b31a35ad7bf0)

 or 
 
* recent release of the Pfam-A database (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/).

To work properly, viralVerify require Prodigal and hmmsearch in your PATH environment variable.


### Optional BLAST verification

You can verify your output by BLAST to check if you found novel viruses or plasmids. In this case, you need to have blastn in your $PATH, Biopython installed, and provide a path to the nucleotide database (e.g. local copy of the NCBI nt database). For each contig we report information (e-value, query coverage, identity and subject title) about its best blast hit in the provided database.


### Usage 

    viralverify 
            -f Input fasta file
            -o output_directory 
            --hmm HMM  Path to HMM database

            Optional arguments:
            -h, --help  Show the help message and exit
            --db DB     Run BLAST on input contigs against provided db
            -t          Number of threads
            -thr THR    Sensitivity threshold (minimal absolute score to classify sequence, default = 7)
            -p          Output predicted plasmidic contigs separately


Output file: comma-separated table *<input_file>_result_table.csv*

Output format: contig name, prediction result, log-likelihood ratio, list of predicted HMMs
  
Fasta files with prediction results can be found in the *Prediction_results_fasta* folder
  
To decrease number of false positives (at the expense of potential false negatives) you may increase the detection threshold, provided as an optional argument.

### Retraining classifier

You can retrain the classifier with your custom data using provided *training_script*. It takes viral, chromosomal and plasmid (optionally) training sequences in fasta format and set of HMMS, predict genes and HMM hits, and returns the frequency table. To use the retrained classifier, replace the "classifier_table.txt" file in the viralVerify directory with the obtained table.
