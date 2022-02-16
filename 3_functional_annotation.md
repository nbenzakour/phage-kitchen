## Phage functional annotation

### BALROG - Bacterial Annotation by Learned Representation Of Genes

**Description:**

https://github.com/salzberg-lab/Balrog
https://www.biorxiv.org/content/10.1101/2020.09.06.285304v1

Balrog is a prokaryotic gene finder based on a Temporal Convolutional Network. We took a data-driven approach to prokaryotic gene finding, relying on the large and diverse collection of already-sequenced genomes. By training a single, universal model of bacterial genes on protein sequences from many different species, we were able to match the sensitivity of current gene finders while reducing the overall number of gene predictions. Balrog does not need to be refit on any new genome.


### DRAM - Distilled and Refined Annotation of Metabolism

**Description:**

https://github.com/shafferm/DRAM

DRAM (Distilled and Refined Annotation of Metabolism) is a tool for annotating metagenomic assembled genomes and VirSorter identified viral contigs. DRAM annotates MAGs and viral contigs using KEGG (if provided by the user), UniRef90, PFAM, dbCAN, RefSeq viral, VOGDB and the MEROPS peptidase database as well as custom user databases. DRAM is run in two stages. First an annotation step to assign database identifiers to gene and then a distill step to curate these annotations into useful functional categories. Additionally viral contigs are further analyzed during to identify potential AMGs. This is done via assigning an auxiliary score and flags representing the confidence that a gene is both metabolic and viral.

For more detail on DRAM and how DRAM works please see our paper as well as the wiki.
https://academic.oup.com/nar/article/48/16/8883/5884738
https://github.com/shafferm/DRAM/wiki


### MultiPHATE2 - bioinformatics pipeline for functional annotation of phage isolates

**Description:**

https://github.com/carolzhou/multiPhATE2
MULTIPHATE2 - > http://dx.doi.org/10.1093/g3journal/jkab074
MULTIPHATE - > https://doi.org/10.1093/bioinformatics/btz258
PHANNOTATE -> https://doi.org/10.1093/bioinformatics/btz265

**ABOUT THE MULTI-PHATE PIPELINE DRIVER**
MultiPhATE is a command-line program that runs gene finding and the PhATE annotation code over user-specified phage genomes, then performs gene-by-gene comparisons among the genomes. The multiPhate.py code takes a single argument consisting of a configuration file (hereafter referred to as, multiPhate.config; use the file sample.multiPhate.config as starting point) and uses it to specify annotation parameters. Then, multiPhate.py invokes the PhATE pipeline for each genome. See below for the types of annotations that PhATE performs. If two or more genomes are specified by the user, then multiPhATE will run the CompareGeneProfiles code to identify corresponding genes among the genomes.

**ABOUT THE PHATE PIPELINE**
PhATE is a fully automated computational pipeline for identifying and annotating phage genes in genome sequence. PhATE is written in Python 3.7, and runs on Linux and Mac operating systems. Code execution is controled by a configuration file, which can be tailored to run specific gene finders and to blast sequences against specific phage- and virus-centric data sets, in addition to more generic (genome, protein) data sets. See below for the specific databases that are accommodated. PhATE runs at least one gene finding algorithm, then annotates the genome, gene, and protein sequences using nucleotide and protein blast flavors and a set of fasta sequence databases, and uses hmm searches (phmmer, jackhmmer) against these same fasta databases. It also runs hmmscan against the pVOG and VOG hmm profile databases. If more than one gene finder is run, PhATE will provide a side-by-side comparison of the genes called by each gene caller. The user specifies the preferred gene caller, and the genes and proteins predicted by that caller are annotated using blast against the supporting databases (or, the user may specify one of the comparison gene sets: superset, consensus, or commoncore, for functional annotation). Classification of each protein sequence into a pVOG or VOG group is followed by generation of an alignment-ready fasta file. By convention, genome sequence files end with extension, ".fasta"; gene nucleotide fasta files end with, ".fnt", and cds amino-acid fasta files end with, ".faa".


### Phage Commander, an Application for Rapid Gene Identification in Bacteriophage Genomes Using Multiple Programs

**Description:**

https://github.com/sarah-harris/PhageCommander
https://www.liebertpub.com/doi/full/10.1089/phage.2020.0044

We present the Phage Commander application for rapid identification of bacteriophage genes using multiple gene identification programs. Phage Commander runs a bacteriophage genome sequence through nine gene identification programs (and an additional program for identification of tRNAs) and integrates the results within a single output table. Phage Commander also generates formatted output files for direct export to National Center for Biotechnology Information GenBank or genome visualization programs such as DNA Master. 


### PhagePromoter - Predicting promoters in phage genomes

**Description:**

https://academic.oup.com/bioinformatics/article/35/24/5301/5540317

The growing interest in phages as antibacterial agents has led to an increase in the number of sequenced phage genomes, increasing the need for intuitive bioinformatics tools for performing genome annotation. The identification of phage promoters is indeed the most difficult step of this process. Due to the lack of online tools for phage promoter prediction, we developed PhagePromoter, a tool for locating promoters in phage genomes, using machine learning methods. This is the first online tool for predicting promoters that uses phage promoter data and the first to identify both host and phage promoters with different motifs.

Availability and implementation
This tool was integrated in the Galaxy framework and it is available online at: https://bit.ly/2Dfebfv.


### PhageTerm: a tool for fast and accurate determination of phage termini and packaging mechanism using next-generation sequencing data

**Description:**

https://www.nature.com/articles/s41598-017-07910-5
https://sourceforge.net/projects/phageterm

 In this work, we demonstrate how it is possible to recover more information from sequencing data than just the phage genome. We developed a theoretical and statistical framework to determine DNA termini and phage packaging mechanisms using NGS data. Our method relies on the detection of biases in the number of reads, which are observable at natural DNA termini compared with the rest of the phage genome. We implemented our method with the creation of the software PhageTerm and validated it using a set of phages with well-established packaging mechanisms representative of the termini diversity, i.e. 5′cos (Lambda), 3′cos (HK97), pac (P1), headful without a pac site (T4), DTR (T7) and host fragment (Mu). In addition, we determined the termini of nine Clostridium difficile phages and six phages whose sequences were retrieved from the Sequence Read Archive. PhageTerm is freely available (https://sourceforge.net/projects/phageterm), as a Galaxy ToolShed and on a Galaxy-based server (https://galaxy.pasteur.fr).


### PhageTermVirome - High-throughput identification of viral termini and packaging mechanisms in virome datasets

**Description:**

https://www.nature.com/articles/s41598-021-97867-3

Here, we introduce PhageTermVirome (PTV) as a tool for the easy and rapid high-throughput determination of phage termini and packaging mechanisms using modern large-scale metagenomics datasets. We successfully tested the PTV algorithm on a mock virome dataset and then used it on two real virome datasets to achieve the rapid identification of more than 100 phage termini and packaging mechanisms, with just a few hours of computing time. Because PTV allows the identification of free fully formed viral particles (by recognition of termini present only in encapsidated DNA), it can also complement other virus identification softwares to predict the true viral origin of contigs in viral metagenomics datasets. PTV is a novel and unique tool for high-throughput characterization of phage genomes, including phage termini identification and characterization of genome packaging mechanisms. This software should help researchers better visualize, map and study the virosphere. PTV is freely available for downloading and installation at https://gitlab.pasteur.fr/vlegrand/ptv.


### Phamerator & BYU-Phamerator

**Description:**

**Phamerator**

https://phamerator.org/
https://phagesdb.org/Phamerator/faq/
https://github.com/scresawn/Phamerator
http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=21991981

Phamerator is a comparative genomics and genome exploration tool designed and written by Dr. Steve Cresawn of James Madison University.

In 2017, Phamerator transitioned from a Linux-based program to be a cross-platform web-based program. It is available at https://phamerator.org/

-----------
**BYU-Phamerator**

https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3018-2

Here we describe modifications to the phage comparative genomics software program, Phamerator, provide public access to the code, and include instructions for creating custom Phamerator databases. We further report genomic analysis techniques to determine phage packaging strategies and identification of the physical ends of phage genomes.

Results
The original Phamerator code can be successfully modified and custom databases can be generated using the instructions we provide. Results of genome map comparisons within a custom database reveal obstacles in performing the comparisons if a published genome has an incorrect complementarity or an incorrect location of the first base of the genome, which are common issues in GenBank-downloaded sequence files. To address these issues, we review phage packaging strategies and provide results that demonstrate identification of the genome start location and orientation using raw sequencing data and software programs such as PAUSE and Consed to establish the location of the physical ends of the genome. These results include determination of exact direct terminal repeats (DTRs) or cohesive ends, or whether phages may use a headful packaging strategy. Phylogenetic analysis using ClustalO and phamily circles in Phamerator demonstrate that the large terminase gene can be used to identify the phage packaging strategy and thereby aide in identifying the physical ends of the genome.


### PhANNs - a fast and accurate tool and web server to classify phage structural proteins

**Description:**

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007845
http://edwards.sdsu.edu/phanns
https://github.com/Adrian-Cantu/PhANNs

PhANNs is a tool to classify any phage ORF as one of 10 structural protein class, or as "others". It uses an ensemble of Artificial Neural Networks. PhANNs predicts the structural class of a phage ORF by running an artificial neural network ensemble that we created against a fasta file of protein sequences. If you upload a multi-fasta file, we’ll provide you estimates of the structural classes of all the proteins, and we’ll let you download the sequences for each class as a fasta file.

---------------------

For any given bacteriophage genome or phage-derived sequences in metagenomic data sets, we are unable to assign a function to 50–90% of genes, or more. Structural protein-encoding genes constitute a large fraction of the average phage genome and are among the most divergent and difficult-to-identify genes using homology-based methods. To understand the functions encoded by phages, their contributions to their environments, and to help gauge their utility as potential phage therapy agents, we have developed a new approach to classify phage ORFs into ten major classes of structural proteins or into an “other” category. The resulting tool is named PhANNs (Phage Artificial Neural Networks). We built a database of 538,213 manually curated phage protein sequences that we split into eleven subsets (10 for cross-validation, one for testing) using a novel clustering method that ensures there are no homologous proteins between sets yet maintains the maximum sequence diversity for training. An Artificial Neural Network ensemble trained on features extracted from those sets reached a test F1-score of 0.875 and test accuracy of 86.2%. PhANNs can rapidly classify proteins into one of the ten structural classes or, if not predicted to fall in one of the ten classes, as “other,” providing a new approach for functional annotation of phage proteins. PhANNs is open source and can be run from our web server or installed locally.


### PHANOTATE

**Description:**

https://academic.oup.com/bioinformatics/article/35/22/4537/5480131
https://github.com/deprekate/PHANOTATE

PHANOTATE is a tool to annotate phage genomes. It uses the assumption that non-coding bases in a phage genome is disadvantageous, and then populates a weighted graph to find the optimal path through the six frames of the DNA where open reading frames are beneficial paths, while gaps and overlaps are penalized paths.


### PHASTER

**Description:**

https://phaster.ca/
http://www.ncbi.nlm.nih.gov/pubmed/27141966

PHASTER (PHAge Search Tool Enhanced Release) is a significant upgrade to the popular PHAST web server for the rapid identification and annotation of prophage sequences within bacterial genomes and plasmids. While the steps in the phage identification pipeline in PHASTER remain largely the same as in the original PHAST, numerous software improvements and significant hardware enhancements have now made PHASTER faster, more efficient, more visually appealing and much more user friendly. In particular, PHASTER is now 4.3X faster than PHAST when analyzing a typical bacterial genome. More specifically, software optimizations have made the backend of PHASTER 2.7X faster than PHAST. Likewise, the addition of more than 120 CPUs to the PHASTER compute cluster have greatly reduced processing times. PHASTER can now process a typical bacterial genome in 3 minutes from the raw sequence alone, or in 1.5 minutes when given a pre-annotated GenBank file. A number of other optimizations have been implemented, including automated algorithms to reduce the size and redundancy of PHASTER’s databases, improvements in handling multiple (metagenomic) queries and high user traffic, and the ability to perform automated look-ups against >14,000 previously PHAST/PHASTER annotated bacterial genomes (which can lead to complete phage annotations in seconds as opposed to minutes). PHASTER’s web interface has also been entirely rewritten. A new graphical genome browser has been added, gene/genome visualization tools have been improved, and the graphical interface is now more modern, robust, and user-friendly.


### PHROKKA

**Description:**

https://github.com/gbouras13/phrokka

phrokka is a fast phage annotation pipeline.

phrokka uses Phanotate (McNair et al 2019 doi:10.1093/bioinformatics/btz265) to conduct gene calling and tRNAscan-SE 2 (Chan et al 2021 https://doi.org/10.1093/nar/gkab688) to call tRNAs.

phrokka then uses the lightweight PHROGS database (https://phrogs.lmge.uca.fr Terzian et al 2021 https://doi.org/10.1093/nargab/lqab067) to conduct annotation. Specifically, each gene is compared against the entire PHROGS database using mmseqs2.

---
phrokka creates a number of output files in different formats.

The 2 main files phrokka generates is phrokka.gff, which is a gff3 format file including the fasta following the gff table annotations.

phrokka also creates phrokka.tbl, which is a flat-file table suitable to be unploaded to the NCBI's Bankit.


### Prophage Hunter - an integrative hunting tool for active prophages

**Description:**

https://academic.oup.com/nar/article/47/W1/W74/5494712

We present Prophage Hunter, a tool aimed at hunting for active prophages from whole genome assembly of bacteria. Combining sequence similarity-based matching and genetic features-based machine learning classification, we developed a novel scoring system that exhibits higher accuracy than current tools in predicting active prophages on the validation datasets. The option of skipping similarity matching is also available so that there's higher chance for novel phages to be discovered. Prophage Hunter provides a one-stop web service to extract prophage genomes from bacterial genomes, evaluate the activity of the prophages, identify phylogenetically related phages, and annotate the function of phage proteins. Prophage Hunter is freely available at https://pro-hunter.bgi.com/.


### RAST

**Description:**

https://link.springer.com/protocol/10.1007/978-1-4939-7343-9_17
https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3965101/


### STEP3

**Description:**

https://step3.erc.monash.edu/
https://journals.asm.org/doi/10.1128/mSystems.00242-21

**A machine-learning approach to define the component parts of bacteriophage virions**
Bacteriophages (phages) are currently under consideration as a means to treat a wide range of bacterial infections, including those caused by drug-resistant “superbugs”. Successful phage therapy protocols require diverse phage in a phage cocktail, with the prospective need to recognize features of diverse phage from under-sampled environments. The effective use of these viruse for therapy depends on a number of factors, not least of which is the sequence-based choices that must be made to identify new phages for development into phage therapy.

Phage virions, i.e. the physical form of the phage that would be delivered to the site of infection, conform to a blue-print that consists of a protein capsid housing the viral genome, and a multicomponent tail. We view these virions as molecular machines, and the machinery of the tail machinery is complex. First and foremost, elements within the tail function to engage a species-specific component on the surface of the host bacterium, thereby initiating the infection cascade. The tail machinery is also responsible for penetrating through the bacterial cell wall, in order that the tip of the tail can enter the bacterial cytoplasm. Then, and only then, is a signal transmitted to the portal at the proximal end of the tail, enabling release of the phage DNA into the tail lumen to permit DNA translocation into the bacterial cell cytoplasm, resulting in bacterial death.

We have developed an ensemble predictor called STEP3 that uses machine-learning algorithms to characterize the components of the machinery in phage virions. STEP3 can be used to understand the universal features of the machinery in phage tails, by accurately classifying proteins with conserved features together into groupings that are not dependent on the ill-considered annotations that currently confuse phage genome data. In the development of STEP3, various types of evolutionary features were sampled, features that were extracted from Position-Specific Scoring Matrix (PSSM), to draw on relationships underpinning the evolutionary history of the various proteins making up the phage virions. Considering the high evolution rates of phage proteins, these features are particularly suitable to detect virion proteins with only distantly related homologies. STEP3 integrated these features into an ensemble framework to achieve a stable and robust prediction performance. The final ensemble model showed a significant improvement in terms of prediction accuracy over current state-of-the-art phage virion protein predictors on extensive 5-fold cross-validation and independent tests.


### VIGA - De novo Viral Genome Annotator

**Description:**

https://github.com/EGTortuero/viga/tree/developer

VIGA is a script written in Python 3 that annotates viral genomes automatically (using a de novo algorithm) and predict the function of their proteins using BLAST and HMMER. This script works in UNIX-based OS, including MacOSX and the Windows Subsystem for Linux.

Programs:

* LASTZ (Harris 2007): it is used to predict the circularity of the contigs. The program is publicly available at https://github.com/lastz/lastz under the MIT licence.
* INFERNAL (Nawrocki and Eddy 2013): it is used to predict ribosomal RNA in the contigs when using the RFAM database (Nawrocki et al. 2015). This program is publicly available at http://eddylab.org/infernal/ under the BSD licence and RFAM database is available at ftp://ftp.ebi.ac.uk/pub/databases/Rfam/
* ARAGORN (Laslett and Canback 2004): it is used to predict tRNA sequences in the contig. This program is publicly available at http://mbio-serv2.mbioekol.lu.se/ARAGORN/ under the GPLv2 licence.
* PILERCR (Edgar 2007): it is used to predict CRISPR repeats in your contig. This program is freely available at http://drive5.com/pilercr/ under a public licence.
* Prodigal (Hyatt et al. 2010): it is used to predict the ORFs. When the contig is smaller than 100,000 bp, MetaProdigal (Hyatt et al. 2012) is automatically activated instead of normal Prodigal. This program is publicly available at https://github.com/hyattpd/prodigal/releases/ under the GPLv3 licence.
* DIAMOND (Buchfink et al. 2015): it is used to predict the function of proteins according to homology. This program is publicly available at https://github.com/bbuchfink/diamond under the GPLv3 licence. Databases must be created from FASTA files according to their instructions before running.
* BLAST+ (Camacho et al. 2008): it is used to predict the function of the predicted proteins according to homology when DIAMOND is not able to retrieve any hit or such hit is a 'hypothetical protein'. This suite is publicly available at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ under the GPLv2 licence. Databases are available at ftp://ftp.ncbi.nlm.nih.gov/blast/db/ or created using makeblastdb command.
* HMMER (Finn et al. 2011): it is used to add more information of the predicted proteins according to Hidden Markov Models. This suite is publicly available at http://hmmer.org/ under the GPLv3 licence. Databases must be in HMM format and an example of potential database is PVOGs (http://dmk-brain.ecn.uiowa.edu/VOG/downloads/All/AllvogHMMprofiles.tar.gz).


### VIRALPRO

**Description:**

http://scratch.proteomics.ics.uci.edu/explanation.html#VIRALpro

VIRALpro is a predictor capable of identifying capsid and tail protein sequences using support vector machines (SVM) with an accuracy estimated to be between 90% and 97%. Predictions are based on the protein amino acid composition, on the protein predicted secondary structure, as predicted by SSpro, and on a boosted linear combination of HMM e-values obtained from 3,380 HMMs built from multiple sequence alignments of specific fragments - called contact fragments - of both capsid and tail sequences.

http://download.igb.uci.edu/