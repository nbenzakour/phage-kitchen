## Phage full pipeline

### Cenote-Taker2: Discover and Annotate Divergent Viral Contigs

**Description:**

https://academic.oup.com/ve/article/7/1/veaa100/6055568
https://github.com/mtisza1/Cenote-Taker2

Cenote-Taker 2 is a dual function bioinformatics tool. On the one hand, Cenote-Taker 2 discovers/predicts virus sequences from any kind of genome or metagenomic assembly. Second, virus sequences/genomes are annotated with a variety of sequences features, genes, and taxonomy. Either the discovery or the the annotation module can be used independently.

Cenote-Taker 2 democratizes virus discovery and sequence annotation.


### CPT Galaxy

**Description:**

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008214

At the Center for Phage Technology (CPT), we developed a suite of phage-oriented tools housed in open, user-friendly web-based interfaces. A Galaxy platform conducts computationally intensive analyses and Apollo, a collaborative genome annotation editor, visualizes the results of these analyses. The collection includes open source applications such as the BLAST+ suite, InterProScan, and several gene callers, as well as unique tools developed at the CPT that allow maximum user flexibility. We describe in detail programs for finding Shine-Dalgarno sequences, resources used for confident identification of lysis genes such as spanins, and methods used for identifying interrupted genes that contain frameshifts or introns. At the CPT, genome annotation is separated into two robust segments that are facilitated through the automated execution of many tools chained together in an operation called a workflow. First, the structural annotation workflow results in gene and other feature calls. This is followed by a functional annotation workflow that combines sequence comparisons and conserved domain searching, which is contextualized to allow integrated evidence assessment in functional prediction. Finally, we describe a workflow used for comparative genomics. Using this multi-purpose platform enables researchers to easily and accurately annotate an entire phage genome. 

The portal can be accessed at https://cpt.tamu.edu/galaxy-pub with accompanying user training material.
https://cpt.tamu.edu/training-material/


### MGnify: the microbiome analysis resource in 2020

**Description:**

https://academic.oup.com/nar/article/48/D1/D570/5614179

MGnify (http://www.ebi.ac.uk/metagenomics) provides a free to use platform for the assembly, analysis and archiving of microbiome data derived from sequencing microbial populations that are present in particular environments. Over the past 2 years, MGnify (formerly EBI Metagenomics) has more than doubled the number of publicly available analysed datasets held within the resource. Recently, an updated approach to data analysis has been unveiled (version 5.0), replacing the previous single pipeline with multiple analysis pipelines that are tailored according to the input data, and that are formally described using the Common Workflow Language, enabling greater provenance, reusability, and reproducibility. MGnify's new analysis pipelines offer additional approaches for taxonomic assertions based on ribosomal internal transcribed spacer regions (ITS1/2) and expanded protein functional annotations. Biochemical pathways and systems predictions have also been added for assembled contigs. MGnify's growing focus on the assembly of metagenomic data has also seen the number of datasets it has assembled and analysed increase six-fold. The non-redundant protein database constructed from the proteins encoded by these assemblies now exceeds 1 billion sequences. Meanwhile, a newly developed contig viewer provides fine-grained visualisation of the assembled contigs and their enriched annotations.


### Phigaro: high throughput prophage sequence annotation

**Description:**

https://www.biorxiv.org/content/10.1101/598243v1

Summary Phigaro is a standalone command-line application that is able to detect prophage regions taking raw genome and metagenome assemblies as an input. It also produces dynamic annotated “prophage genome maps” and marks possible transposon insertion spots inside prophages. It provides putative taxonomic annotations that can distinguish tailed from non-tailed phages. It is applicable for mining prophage regions from large metagenomic datasets.

Availability Source code for Phigaro is freely available for download at https://github.com/bobeobibo/phigaro along with test data. The code is written in Python.


### VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences

**Description:**

https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0
https://github.com/AnantharamanLab/VIBRANT

*Background*
Viruses are central to microbial community structure in all environments. The ability to generate large metagenomic assemblies of mixed microbial and viral sequences provides the opportunity to tease apart complex microbiome dynamics, but these analyses are currently limited by the tools available for analyses of viral genomes and assessing their metabolic impacts on microbiomes.

*Design*
Here we present VIBRANT, the first method to utilize a hybrid machine learning and protein similarity approach that is not reliant on sequence features for automated recovery and annotation of viruses, determination of genome quality and completeness, and characterization of viral community function from metagenomic assemblies. VIBRANT uses neural networks of protein signatures and a newly developed v-score metric that circumvents traditional boundaries to maximize identification of lytic viral genomes and integrated proviruses, including highly diverse viruses. VIBRANT highlights viral auxiliary metabolic genes and metabolic pathways, thereby serving as a user-friendly platform for evaluating viral community function. VIBRANT was trained and validated on reference virus datasets as well as microbiome and virome data.

*Results*
VIBRANT showed superior performance in recovering higher quality viruses and concurrently reduced the false identification of non-viral genome fragments in comparison to other virus identification programs, specifically VirSorter, VirFinder, and MARVEL. When applied to 120,834 metagenome-derived viral sequences representing several human and natural environments, VIBRANT recovered an average of 94% of the viruses, whereas VirFinder, VirSorter, and MARVEL achieved less powerful performance, averaging 48%, 87%, and 71%, respectively. Similarly, VIBRANT identified more total viral sequence and proteins when applied to real metagenomes. When compared to PHASTER, Prophage Hunter, and VirSorter for the ability to extract integrated provirus regions from host scaffolds, VIBRANT performed comparably and even identified proviruses that the other programs did not. To demonstrate applications of VIBRANT, we studied viromes associated with Crohn’s disease to show that specific viral groups, namely Enterobacteriales-like viruses, as well as putative dysbiosis associated viral proteins are more abundant compared to healthy individuals, providing a possible viral link to maintenance of diseased states.

*Conclusions*
The ability to accurately recover viruses and explore viral impacts on microbial community metabolism will greatly advance our understanding of microbiomes, host-microbe interactions, and ecosystem dynamics.

***Completeness estimation**
Scaffold completeness is determined based on four metrics: circularization of scaffold sequence, VOG annotations, total VOG nucleotide replication proteins, and total VOG viral hallmark proteins (Additional File 13: Table S13). In order to be considered a complete genome, a sequence must be identified as likely circular. A kmer-based approach is used to do this. Specifically, the first 20 nucleotides are compared to 20-mer sliding windows within the last 900 bp of the sequence. If a complete match is identified, the sequence is considered a circular template. Scaffolds can also be considered a low-, medium-, or high-quality draft. To benchmark completeness, 2466 NCBI RefSeq viruses identified as Caudovirales, limited to 10 kb in length, were used to estimate completeness by stepwise removing 10% viral sequence at a time. VIBRANT was found to identify 2465 of the 2466 viruses. This set of viruses was additionally used to assess the error rate of cutting provirus regions. Viral genome diagrams to depict genome quality and completeness, provirus predictions, and novel virus identification were made using Geneious Prime 2019.0.3.

### VIRify

**Description:**

https://github.com/EBI-Metagenomics/emg-viral-pipeline

VIRify is a recently developed pipeline for the detection, annotation, and taxonomic classification of viral contigs in metagenomic and metatranscriptomic assemblies. The pipeline is part of the repertoire of analysis services offered by MGnify. VIRify’s taxonomic classification relies on the detection of taxon-specific profile hidden Markov models (HMMs), built upon a set of 22,014 orthologous protein domains and referred to as ViPhOGs.

The pipeline is implemented and available in CWL and Nextflow.


### What the Phage: A scalable workflow for the identification and analysis of phage sequences

**Description:**

https://github.com/replikation/What_the_Phage
https://doi.org/10.1101/2020.07.24.219899

* WtP is a scalable and easy-to-use workflow for phage identification and analysis. Our tool currently combines 10 established phage identification tools
* An attempt to streamline the usage of various phage identification and prediction tools
* The main focus is stability and data filtering/analysis for the user
* The tool is intended for fasta and fastq reads to identify phages in contigs/reads
* Proper prophage detection is not implemented (yet) - but a handful of tools report them - so they are mostly identified