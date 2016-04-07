        =======  QuasiPycies =======

Python 3 module to calculate diversity indices from NGS data. It consists of two main classes (Clusters and Pileup) and extra functions for convertion and parsing different file formats.

class Clusters 
	
	Functions:  shannon_entropy
				simpsons_index
				simpsons_index_of_diversity
				simpsons_Reciprocal_Index
				clusters_object_to_fasta
				filtering_clusters_by_percentage
				number_haplotypes
				

class Pileup
			
	Functions:  nucleotide_diversity
				polymorphic_sites


Extra functions
				filtering_clusters_by_percentage
				sam_to_fasta
				fastq_to_fasta
				multiline_fasta_to_singleline
				clusters_to_single_reads


             ======= Generating input files ========

Input file for Clusters Class
-----------------------------
For the Clusters class, you need to provide a fasta file with sequences representing the clusters. The fasta file must have number of reads	("size=") on each sequence name.

eg:

	>ReadA;size=180
	AAATTTCCCGGG
	>ReadB;size=20
	TTTAAATTTGGG 

The easieast way to do it is using Usearch (http://www.drive5.com/usearch/download.html) and the command cluster_fast:

$ ./usearch -cluster_fast [fasta input] -id 1.0 -centroids [fasta output] -sizeout

The id option specifies the percentage of identity needed for clustering the sequences. I use 1.0 which means that clusters are formed only by identical sequences.

The sizeout option prints the number of reads for each cluster on the end of the sequence name (;size=).

Next step is to use the function multiline_fasta_to_singleline from the quasispycies module. The reason for it is because the output fasta from usearch cluster_fast is fasta file with the nucleotide sequence information in multiple lines instead in just one line. Most of the scripts on the quasispycies were designed to work with single line fasta, otherwise it would just use the informaton from the first line of each nucleotide sequence.

Finally I recommend filtering the clusters by percentage based on the error rate of your sequencing and PCR. You can do that by sequencing a viral stock or plasmid and checking the percentage of variants. In my experience, using clusters above 1% is a conservative value when you are dealing with Illumina MiSeq amplicon sequencing. To do it you can use the function filtering_clusters_by_percentage. To use it you have to provide the name of the fasta file (single line fasta), the percentage to filter (if you want to filter only clusters above 1% you should write 0.01) and the fasta output file name. 


Input file for Pileup Class
---------------------------
The input file for the Pileup Class must be in mpileup format. The mpileup file can be generated from a sam/bam file using samtools. The information of stacked bases for each position in the alignment is used to calculated intrapopulation 	heterogeneity. To do it first you need to create a sam file using your reads mapped to a reference. You can use any mapping tool, but I would advise you to use bbmap (https://sourceforge.net/projects/bbmap/), because it is a fast and consistent sequence mapper.

$. bbmap/bbmap.sh ref=[reference.fasta] in1=[input.fastq] out=[output.sam] outputunmapped=f minid=0.6 qtrim=t trimq=20 maq=20 ordered=t -Xmx3g

It is important to use a reference sequence with trimmed primer regions. Otherwise when you generate a mpileup from your alignment, you would get no information from the bases on the primer region because you trimmed the reads before.

The output from bbmap must be a sam file. Next step is to use Samtools to generate the mpileup file. First you need to convert the sam file to bam.

$ samtools view -S -b [file.sam] > [file.bam]

And sort bam file.

For Samtools 1.3:
$ samtools sort -o [sorted.bam] [file.bam]

For Samtools 1.2:
$ samtools sort [file.bam] [file_sorted]

Next, you need to index the sorted bam file:

$ samtools index -b [sorted.bam]

And finally generate the mpileup file:
samtools mpileup -B -d 10000 -f [reference.fasta] [sorted_indexed.bam] -o [output.mpileup]

For more information about Samtools please access (http://www.htslib.org/doc/samtools.html)


Using the quasispycies module
-----------------------------
I recommend you to use a Python notebook system like iPython (https://ipython.org/) to work with your data. IPython is much more flexible than the regular Python IDLE and it helps you store, organize and test your blocks of code in a much more consistent way.

To use the quasispycies module you can save it in your working directory or copy it to the folder where the other packages and modules are stored on your computer. To get the path you can type on your Python3 environment:

$ import sys
$ sys.path

To load the quasispycies module you can type:

$ import quasispycies as qp

And you would have to use qp. preceding each function or class you use.

eg.

$ pileup_A = qp.Pileup("mpileup_file_A")
$ qp.fastq_to_fasta("fastq_file", "fasta_file")

Otherwise you can use:

$ from quasispycies import *

And all the classes and functions would be loaded from the module and it would not be necessary to use the qp. preceding the calls.


Brief description of the classes and functions
----------------------------------------------

To get more details about the classes and functions in the module, you can type help() and the respective name inside the brackets to print the docstrings.

eg.

$ help(Clusters)
$ help(sam_to_fasta)


























