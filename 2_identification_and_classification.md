## Phage identification & classification

### Cenote_Unlimited_Breadsticks

**Description:**

https://github.com/mtisza1/Cenote_Unlimited_Breadsticks

Unlimited Breadsticks uses probabilistic models (i.e. HMMs) of virus hallmark genes to identify virus sequences from any dataset of contigs (e.g. metagenomic assemblies) or genomes (e.g. bacterial genomes). Optionally, Unlimited Breadsticks will use gene content information to remove flanking cellular chromosomes from contigs representing putative prophages. Generally, the prophage-cellular chromosome boundary will be identified within 100 nt - 2000 nt of the actual location.

+ The code is currently functional. Feel free to consume Unlimited Breadsticks at will.
+ Minor update to handle very large contig files AND update to HMM databases on June 16th, 2021

Unlimited Breadsticks is derived from Cenote-Taker 2, but several time-consuming computations are skipped in order to analyze datasets as quickly as possible. Also, Unlimited Breadsticks only takes approximately 16 minutes to download and install (Cenote-Taker 2 takes about 2 hours due to large databases required for thorough sequence annotation). See installation instructions below.

**Limitations**

Compared to Cenote-Taker 2, there are a few limitations.

* Unlimited Breadsticks does not do post-hallmark-gene-identification computations to flag plasmid and conjugative element sequences that occasionally slip through.
* Unlimited Breadsticks does not make genome maps for manual inspection of putative viruses.
* Contigs are not extensively annotated by Unlimited Breadsticks. No genome maps are created.


### CheckV - assesses the quality and completeness of metagenome-assembled viral genomes

**Description:**

https://www.nature.com/articles/s41587-020-00774-7
https://portal.nersc.gov/CheckV/
https://bitbucket.org/berkeleylab/CheckV

Here we present CheckV, an automated pipeline for identifying closed viral genomes, estimating the completeness of genome fragments and removing flanking host regions from integrated proviruses. CheckV estimates completeness by comparing sequences with a large database of complete viral genomes, including 76,262 identified from a systematic search of publicly available metagenomes, metatranscriptomes and metaviromes. After validation on mock datasets and comparison to existing methods, we applied CheckV to large and diverse collections of metagenome-assembled viral sequences, including IMG/VR and the Global Ocean Virome.


### DeePhage: a tool for identifying temperate phage-derived and virulent phage-derived sequence in metavirome data using deep learning

**Description:**

DeePhage is designed to identify metavirome sequences as temperate phage-derived and virulent phage-derived sequences. The program calculate a score reflecting the likelihood of each input fragment as temperate phage-derived and virulent phage-derived sequences. DeePhage can run either on the virtual machine or physical host. For non-computer professionals, we recommend running the virtual machine version of DeePhage on local PC. In this way, users do not need to install any dependency package. If GPU is available, you can also choose to run the physical host version. This version can automatically speed up with GPU and is more suitable to handle large scale data. The program is also available at http://cqb.pku.edu.cn/ZhuLab/DeePhage/.


### DeepVirFinder: Identifying viruses from metagenomic data by deep learning

**Description:**

https://github.com/jessieren/DeepVirFinder
https://link.springer.com/article/10.1007/s40484-019-0187-4
https://link.springer.com/content/pdf/10.1007/s40484-019-0187-4.pdf

DeepVirFinder predicts viral sequences using deep learning method. The method has good prediction accuracy for short viral sequences, so it can be used to predict sequences from the metagenomic data.

DeepVirFinder significantly improves the prediction accuracy compared to our k-mer based method VirFinder by using convolutional neural networks (CNN). CNN can automatically learn genomic patterns from the viral and prokaryotic sequences and simultaneously build a predictive model based on the learned genomic patterns. The learned patterns are represented in the form of weight matrices of size 4 by k, where k is the length of the pattern. This representation is similar to the position weight matrix (PWM), the commonly used representation of biological motifs, which are also of size 4 by k and each column specifies the probabilities of having the 4 nucleotides at that position. When only one type of nucleotide can be chosen at each position with probability 1, the motif degenerates to a k-mer. Thus, the CNN is a natural generalization of k-mer based model. The more flexible CNN model indeed outperforms the k-mer based model on viral sequence prediction problem.


### [75] Demovir

**Created:** 09/11/2021 15:34:37

**Created by:** nbenzakour



**Description:**

https://github.com/feargalr/Demovir

Democratic taxonomic classification of viral contigs to Order and Family level

When performing metagenomic sequencing of Viral-Like Particle (VLP), the majority of returned sequences often bare little to no homology to reference sequences - Viral Dark Matter. Frequently it may be useful to know which viral taxonomic group these novel viruses are likely to belong to as this will give information about nucleic acid type, size and behaviour.

DemoVir will classify viral contigs to the Order or Family taxonomic level by comparing genes on the amino acid level against the viral subset of the TrEMBL database, and then taking a vote of the Order and Family hits. Homology searches are performed by Usearch in order to increase speed. This type of method has previously been implemented in multiple published virome studies but to our knowledge none have performed benchmarking or made it available as a simple executable script easily downloaded and installed.

**Note - DemoVir is for classification of sequences into viral families and orders only and should not be used for discrminating viral contigs from bacterial/archael/eukaryotic sequences in a metagenomic sample.**

### DEPtH - Detection and Extraction of Phages Tool

**Description:**

DEPhT is a new tool for identifying prophages in bacteria, and was developed with a particular interest in being able to rapidly scan hundreds to thousands of genomes and accurately extract complete (likely active) prophages from them.

A detailed manuscript has been submitted to Nucleic Acids Research, but in brief DEPhT works by using genome architecture (rather than homology) to identify genomic regions likely to contain a prophage. Any regions with phage-like architecture (characterized as regions with high gene density and few transcription direction changes) are then further scrutinized using two passes of homology detection. 

* The first pass identifies genes on putative prophages that are homologs of (species/clade/genus-level) conserved bacterial genes, and uses any such genes to disrupt the prophage prediction. 
* The second pass (disabled in the 'fast' runmode) identifies genes on putative prophages that are homologs of conserved, functionally annotated phage genes. 
* Finally, prophage regions that got through the previous filters are subjected to a BLASTN-based attL/attR detection scheme that gives DEPhT better boundary detection than any tool we are aware of.


### GRAViTy: Genome Relationships Applied to Virus Taxonomy

**Description:**

http://gravity.cvr.gla.ac.uk/
https://github.com/PAiewsakun/GRAViTy

GRAViTy - "Genome Relationships Applied to Virus Taxonomy" is an analysis pipeline that is effective at reproducing the current assignments of viruses at family level as well as inter-family groupings into Orders (Aiewsakun and Simmonds, 2018). It can additionally be used to correctly differentiate assigned viruses from unassigned viruses and classify them into correct taxonomic groups. The method provides a rapid and objective means to explore metagenomic viral diversity and make informed recommendations for classification that are consistent with the current ICTV taxonomic framework. Methods like GRAViTy are increasingly required as the vast diversity of viruses found in metagenomic sequence datasets is explored.

https://www.microbiologyresearch.org/content/journal/jgv/10.1099/jgv.0.001110
https://link.springer.com/article/10.1007%2Fs00705-018-3938-z
https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0422-7


### Hecatomb

**Description:**

https://hecatomb.readthedocs.io/en/latest/
https://github.com/shandley/hecatomb

A hecatomb is a great sacrifice or an extensive loss. Heactomb the software empowers an analyst to make data driven decisions to 'sacrifice' false-positive viral reads from metagenomes to enrich for true-positive viral reads. This process frequently results in a great loss of suspected viral sequences / contigs.


### INPHARED - INfrastructure for a PHAge REference Database: Identification of large-scale biases in the current collection of phage genomes.

**Description:**

https://doi.org/10.1101/2021.05.01.442102
https://github.com/RyanCook94/inphared

inphared.pl (INfrastructure for a PHAge REference Database) is a perl script which downloads and filters phage genomes from Genbank to provide the most complete phage genome database possible.

Providing up-to-date bacteriophage genome databases, metrics and useful input files for a number of bioinformatic pipelines including vConTACT2 and MASH. The aim is to produce a useful starting point for viral genomics and meta-omics.

https://leicester.figshare.com/articles/dataset/INPHARED_DATABASE/14242085


### IMGVR v3: an integrated ecological and evolutionary framework for interrogating genomes of uncultivated viruses

**Description:**

https://doi.org/10.1093/nar/gkaa946

Viruses are integral components of all ecosystems and microbiomes on Earth. Through pervasive infections of their cellular hosts, viruses can reshape microbial community structure and drive global nutrient cycling. Over the past decade, viral sequences identified from genomes and metagenomes have provided an unprecedented view of viral genome diversity in nature. Since 2016, the IMG/VR database has provided access to the largest collection of viral sequences obtained from (meta)genomes. Here, we present the third version of IMG/VR, composed of 18 373 cultivated and 2 314 329 uncultivated viral genomes (UViGs), nearly tripling the total number of sequences compared to the previous version. These clustered into 935 362 viral Operational Taxonomic Units (vOTUs), including 188 930 with two or more members. UViGs in IMG/VR are now reported as single viral contigs, integrated proviruses or genome bins, and are annotated with a new standardized pipeline including genome quality estimation using CheckV, taxonomic classification reflecting the latest ICTV update, and expanded host taxonomy prediction. The new IMG/VR interface enables users to efficiently browse, search, and select UViGs based on genome features and/or sequence similarity. 

IMG/VR v3 is available at https://img.jgi.doe.gov/vr, and the underlying data are available to download at https://genome.jgi.doe.gov/portal/IMG_VR.


### MARVEL, a Tool for Prediction of Bacteriophage Sequences in Metagenomic Bins

**Description:**

https://www.frontiersin.org/articles/10.3389/fgene.2018.00304/full

Here we present MARVEL, a tool for prediction of double-stranded DNA bacteriophage sequences in metagenomic bins. MARVEL uses a random forest machine learning approach. 

https://github.com/LaboratorioBioinformatica/MARVEL


### MetaPhinder—Identifying Bacteriophage Sequences in Metagenomic Data Sets

**Description:**

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163111

Bacteriophages are the most abundant biological entity on the planet, but at the same time do not account for much of the genetic material isolated from most environments due to their small genome sizes. They also show great genetic diversity and mosaic genomes making it challenging to analyze and understand them. Here we present MetaPhinder, a method to identify assembled genomic fragments (i.e.contigs) of phage origin in metagenomic data sets. The method is based on a comparison to a database of whole genome bacteriophage sequences, integrating hits to multiple genomes to accomodate for the mosaic genome structure of many bacteriophages. The method is demonstrated to out-perform both BLAST methods based on single hits and methods based on k-mer comparisons. MetaPhinder is available as a web service at the Center for Genomic Epidemiology https://cge.cbs.dtu.dk/services/MetaPhinder/, while the source code can be downloaded from https://bitbucket.org/genomicepidemiology/metaphinder or https://github.com/vanessajurtz/MetaPhinder.

### METAVIRALSPADES: assembly of viruses from metagenomic data

**Description:**

https://academic.oup.com/bioinformatics/article/36/14/4126/5837667

We describe a METAVIRALSPADES tool for identifying viral genomes in metagenomic assembly graphs that is based on analyzing variations in the coverage depth between viruses and bacterial chromosomes. We benchmarked METAVIRALSPADES on diverse metagenomic datasets, verified our predictions using a set of virus-specific Hidden Markov Models and demonstrated that it improves on the state-of-the-art viral identification pipelines.

**Availability and implementation**
METAVIRALSPADES includes VIRALASSEMBLY, VIRALVERIFY and VIRALCOMPLETE modules that are available as standalone packages:
https://github.com/ablab/spades/tree/metaviral_publication,
https://github.com/ablab/viralVerify/ and
https://github.com/ablab/viralComplete/.


### MGV - Metagenomic Gut Virus catalogue viral detection pipeline

**Description:**

https://github.com/snayfach/MGV
https://www.nature.com/articles/s41564-021-00928-6

The Metagenomic Gut Virus catalogue improves detection of viruses in stool metagenomes and accounts for nearly 40% of CRISPR spacers found in human gut Bacteria and Archaea. We also produced a catalogue of 459,375 viral protein clusters to explore the functional potential of the gut virome

Access to the full catalogue of viral genomes, protein clusters, diversity-generating retroelements and CRISPR spacers is provided without restrictions at https://portal.nersc.gov/MGV. Any requests for further data should be directed to the corresponding authors.


### PhaGCN - GCN based model classifier

**Description:**

https://academic.oup.com/bioinformatics/article/37/Supplement_1/i25/6319660
https://github.com/KennthShang/PhaGCN

PhaGCN is a GCN based model, which can learn the species masking feature via deep learning classifier, for new Phage taxonomy classification. To use PhaGCN, you only need to input your contigs to the program.


### PPR-Meta: a tool for identifying phages and plasmids from metagenomic fragments using deep learning

**Description:**

https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6586199/

We present PPR-Meta, a 3-class classifier that allows simultaneous identification of both phage and plasmid fragments from metagenomic assemblies. PPR-Meta consists of several modules for predicting sequences of different lengths. Using deep learning, a novel network architecture, referred to as the Bi-path Convolutional Neural Network, is designed to improve the performance for short fragments. PPR-Meta demonstrates much better performance than currently available similar tools individually for phage or plasmid identification, while testing on both artificial contigs and real metagenomic data. PPR-Meta is freely available via http://cqb.pku.edu.cn/ZhuLab/PPR_Meta or https://github.com/zhenchengfang/PPR-Meta.


### Seeker: alignment-free identification of bacteriophage genomes by deep learning

**Description:**

https://academic.oup.com/nar/article/48/21/e121/5921300

Recent advances in metagenomic sequencing have enabled discovery of diverse, distinct microbes and viruses. Bacteriophages, the most abundant biological entity on Earth, evolve rapidly, and therefore, detection of unknown bacteriophages in sequence datasets is a challenge. Most of the existing detection methods rely on sequence similarity to known bacteriophage sequences, impeding the identification and characterization of distinct, highly divergent bacteriophage families. Here we present Seeker, a deep-learning tool for alignment-free identification of phage sequences. Seeker allows rapid detection of phages in sequence datasets and differentiation of phage sequences from bacterial ones, even when those phages exhibit little sequence similarity to established phage families. We comprehensively validate Seeker's ability to identify previously unidentified phages, and employ this method to detect unknown phages, some of which are highly divergent from the known phage families. We provide a web portal (seeker.pythonanywhere.com) and a user-friendly Python package (github.com/gussow/seeker) allowing researchers to easily apply Seeker in metagenomic studies, for the detection of diverse unknown bacteriophages.


### vCONTACT2 - Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks

**Description:**

https://bitbucket.org/MAVERICLab/vcontact2/wiki/Home
https://www.nature.com/articles/s41587-019-0100-8

We present vConTACT v.2.0, a network-based application utilizing whole genome gene-sharing profiles for virus taxonomy that integrates distance-based hierarchical clustering and confidence scores for all taxonomic predictions. We report near-identical (96%) replication of existing genus-level viral taxonomy assignments from the International Committee on Taxonomy of Viruses for National Center for Biotechnology Information virus RefSeq. Application of vConTACT v.2.0 to 1,364 previously unclassified viruses deposited in virus RefSeq as reference genomes produced automatic, high-confidence genus assignments for 820 of the 1,364. We applied vConTACT v.2.0 to analyze 15,280 Global Ocean Virome genome fragments and were able to provide taxonomic assignments for 31% of these data, which shows that our algorithm is scalable to very large metagenomic datasets. Our taxonomy tool can be automated and applied to metagenomes from any environment for virus classification.

---

Version 1

https://peerj.com/articles/3243/
vConTACT: an iVirus tool to classify double-stranded DNA viruses that infect Archaea and Bacteria 


### VICTOR: genome-based phylogeny and classification of prokaryotic viruses - online only?

**Description:**

https://doi.org/10.1093/bioinformatics/btx440
https://ggdc.dsmz.de/victor.php

We here present a novel in silico framework for phylogeny and classification of prokaryotic viruses, in line with the principles of phylogenetic systematics, and using a large reference dataset of officially classified viruses. The resulting trees revealed a high agreement with the classification. Except for low resolution at the family level, the majority of taxa was well supported as monophyletic. Clusters obtained with distance thresholds chosen for maximizing taxonomic agreement appeared phylogenetically reasonable, too. Analysis of an expanded dataset, containing >4000 genomes from public databases, revealed a large number of novel species, genera, subfamilies and families.


### ViPTree : the Viral Proteomic Tree server

**Description:**

https://www.genome.jp/viptree/

The ViPTree server generates a "proteomic tree" of viral genome sequences based on genome-wide sequence similarities computed by tBLASTx. The original proteomic tree concept (i.e., "the Phage Proteomic Tree”) was developed by Rohwer and Edwards, 2002. A proteomic tree is a dendrogram that reveals global genomic similarity relationships between tens, hundreds, and thousands of viruses. It has been shown that viral groups identified in a proteomic tree well correspond to established viral taxonomies. The proteomic tree approach is effective to investigate genomes of newly sequenced viruses as well as those identified in metagenomes.

2021-10-04 version 1.9.1
Version of Virus-Host DB: RefSeq release 207

https://github.com/yosuken/ViPTreeGen

ViPTreeGen is a tool for automated generation of viral "proteomic tree" by computing genome-wide sequence similarities based on tBLASTx results.


### VirClust – a tool for hierarchical clustering, core gene detection and annotation of (prokaryotic) viruses

**Description:**

https://doi.org/10.1101/2021.06.14.448304

Here, VirClust is presented – a novel tool capable of performing 

* hierarchical clustering of viruses based on intergenomic distances calculated from their protein cluster content, 
* identification of core proteins and 
* annotation of viral proteins. VirClust groups proteins into clusters both based on BLASTP sequence similarity, which identifies more related proteins, and also based on hidden markow models (HMM), which identifies more distantly related proteins. 

Furthermore, VirClust provides an integrated visualization of the hierarchical clustering tree and of the distribution of the protein content, which allows the identification of the genomic features responsible for the respective clustering. By using different intergenomic distances, the hierarchical trees produced by VirClust can be split into viral genome clusters of different taxonomic ranks. VirClust is freely available, as web-service (virclust.icbm.de) and stand-alone tool.


### VirFinder: a novel k-mer based tool for identifying viral sequences from assembled metagenomic data

**Description:**

https://link.springer.com/epdf/10.1186/s40168-017-0283-5?
https://github.com/jessieren/VirFinder

**Background:** Identifying viral sequences in mixed metagenomes containing both viral and host contigs is a critical first step in analyzing the viral component of samples. Current tools for distinguishing prokaryotic virus and host contigs primarily use gene-based similarity approaches. Such approaches can significantly limit results especially for short contigs that have few predicted proteins or lack proteins with similarity to previously known viruses. 

**Methods:** We have developed VirFinder, the first k-mer frequency based, machine learning method for virus contig identification that entirely avoids gene-based similarity searches. VirFinder instead identifies viral sequences based on our empirical observation that viruses and hosts have discernibly different k-mer signatures. VirFinder’s performance in correctly identifying viral sequences was tested by training its machine learning model on sequences from host and viral genomes sequenced before 1 January 2014 and evaluating on sequences obtained after 1 January 2014. 

**Results:** VirFinder had significantly better rates of identifying true viral contigs (true positive rates (TPRs)) than VirSorter, the current state-of-the-art gene-based virus classification tool, when evaluated with either contigs subsampled from complete genomes or assembled from a simulated human gut metagenome. For example, for contigs subsampled from complete genomes, VirFinder had 78-, 2.4-, and 1.8-fold higher TPRs than VirSorter for 1, 3, and 5 kb contigs, respectively, at the same false positive rates as VirSorter (0, 0.003, and 0.006, respectively), thus VirFinder works considerably better for small contigs than VirSorter. VirFinder furthermore identified several recently sequenced virus genomes (after 1 January 2014) that VirSorter did not and that have no nucleotide similarity to previously sequenced viruses, demonstrating VirFinder’s potential advantage in identifying novel viral sequences. Application of VirFinder to a set of human gut metagenomes from healthy and liver cirrhosis patients reveals higher viral diversity in healthy individuals than cirrhosis patients. We also identified contig bins containing crAssphage-like contigs with higher abundance in healthy patients and a putative Veillonella genus prophage associated with cirrhosis patients. 

**Conclusions:** This innovative k-mer based tool complements gene-based approaches and will significantly improve prokaryotic viral sequence identification, especially for metagenomic-based studies of viral ecology. 

**Keywords:** Metagenome, Virus, k-mer, Human gut, Liver cirrhosis


### VIRIDIC (Virus Intergenomic Distance Calculator) computes pairwise intergenomic distancessimilarities amongst viral genomes.

**Description:**

http://rhea.icbm.uni-oldenburg.de/VIRIDIC/
https://doi.org/10.3390/v12111268

VIRIDIC stand-alone is available now (see download tab). You can use it for jobs with high computational demand and/or for implementing it in your own pipelines. It is very easy to install on your own servers (it is wrapped as a Singularity). You can continue to use the VIRIDIC web-service for small to medium projects (e.g. up to 200 phages per project, no viromes please, they will crash our resources and the analysis will fail).


### VirNET - A deep attention model for viral reads identification

**Description:**

https://github.com/alyosama/virnet
VirNet: A deep attention model for viral reads identification

This tool is able to identifiy viral sequences from a mixture of viral and bacterial sequences. Also, it can purify viral metagenomic data from bacterial contamination


### VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses

**Description:**

https://github.com/jiarong/VirSorter2
https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y

VirSorter2 applies a multi-classifier, expert-guided approach to detect diverse DNA and RNA virus genomes. It has made major updates to its previous version:

* work with more viral groups including dsDNA phages, ssDNA viruses, RNA viruses, NCLDV (Nucleocytoviricota), lavidaviridae (virophages); 
* apply machine learning to estimate viralness using genomic features including structural/functional/taxonomic annotation and viral hallmark genes; 
* train with high quality virus genomes from metagenomes or other sources.

A tutorial/SOP on how to quality control VirSorter2 results is available https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-btv8nn9w

Source code of VirSorter2 is freely available (https://bitbucket.org/MAVERICLab/virsorter2), and VirSorter2 is also available both on bioconda and as an iVirus app on CyVerse (https://de.cyverse.org/de).


### VPF-Class: taxonomic assignment and host prediction of uncultivated viruses based on viral protein families

**Description:**

https://academic.oup.com/bioinformatics/article-abstract/37/13/1805/6104829
https://github.com/biocom-uib/vpf-tools
Supplementary information http://bioinfo.uib.es/~recerca/VPF-Class/

Viral Protein Families (VPFs) can be used for the robust identification of new viral sequences in large metagenomics datasets. Despite the importance of VPF information for viral discovery, VPFs have not yet been explored for determining viral taxonomy and host targets.