## Phage lifestyle & host

### BACPHLIP

**Description:**

https://pubmed.ncbi.nlm.nih.gov/33996289/

Bacteriophages are broadly classified into two distinct lifestyles: temperate and virulent. Temperate phages are capable of a latent phase of infection within a host cell (lysogenic cycle), whereas virulent phages directly replicate and lyse host cells upon infection (lytic cycle). Accurate lifestyle identification is critical for determining the role of individual phage species within ecosystems and their effect on host evolution. Here, we present BACPHLIP, a BACterioPHage LIfestyle Predictor. BACPHLIP detects the presence of a set of conserved protein domains within an input genome and uses this data to predict lifestyle via a Random Forest classifier that was trained on a dataset of 634 phage genomes. On an independent test set of 423 phages, BACPHLIP has an accuracy of 98% greatly exceeding that of the previously existing tools (79%). BACPHLIP is freely available on GitHub 
(https://github.com/adamhockenberry/bacphlip) 
and the code used to build and test the classifier is provided in a separate repository 
(https://github.com/adamhockenberry/bacphlip-model-dev)
 for users wishing to interrogate and re-train the underlying classification model.

### DeepHost

**Description:**

https://github.com/deepomicslab/DeepHost

DeepHost is a phage host prediction tool.


### HTP - Host Taxon Predictor

**Description:**

https://www.nature.com/articles/s41598-019-39847-2

https://github.com/wojciech-galan/Viral_feature_extractor
The initial repo was split into two parts. This part contains a software designed to fetch complete viral genomic reference sequences from NCBI Nucleotide, get viral host's lineage from NCBI Taxonomy and transform the sequence into some features. The second part, available on https://github.com/wojciech-galan/viruses_classifier has been designed to infer host of previously unknown virus.


https://github.com/wojciech-galan/viruses_classifier

Recent advances in metagenomics provided a valuable alternative to culture-based approaches for better sampling viral diversity. However, some of newly identified viruses lack sequence similarity to any of previously sequenced ones, and cannot be easily assigned to their hosts. Here we present a bioinformatic approach to this problem. We developed classifiers capable of distinguishing eukaryotic viruses from the phages achieving almost 95% prediction accuracy. The classifiers are wrapped in Host Taxon Predictor (HTP) software written in Python which is freely available at https://github.com/wojciech-galan/viruses_classifier. HTP’s performance was later demonstrated on a collection of newly identified viral genomes and genome fragments. In summary, HTP is a culture- and alignment-free approach for distinction between phages and eukaryotic viruses. We have also shown that it is possible to further extend our method to go up the evolutionary tree and predict whether a virus can infect narrower taxa.

### MAVRICH - Bacteriophage evolution differs by host, lifestyle and genome

**Description:**

https://www.nature.com/articles/nmicrobiol2017112

Bacteriophages play key roles in microbial evolution, marine nutrient cycling and human disease. Phages are genetically diverse, and their genome architectures are characteristically mosaic, driven by horizontal gene transfer with other phages and host genomes. As a consequence, phage evolution is complex and their genomes are composed of genes with distinct and varied evolutionary histories. However, there are conflicting perspectives on the roles of mosaicism and the extent to which it generates a spectrum of genome diversity or genetically discrete populations. Here, we show that bacteriophages evolve within two general evolutionary modes that differ in the extent of horizontal gene transfer by an order of magnitude. Temperate phages distribute into high and low gene flux modes, whereas lytic phages share only the lower gene flux mode. The evolutionary modes are also a function of the bacterial host and different proportions of temperate and lytic phages are distributed in either mode depending on the host phylum. Groups of genetically related phages fall into either the high or low gene flux modes, suggesting there are genetic as well as ecological drivers of horizontal gene transfer rates. Consequently, genome mosaicism varies depending on the host, lifestyle and genetic constitution of phages.


### PhageAI - PhageAI is an Artificial Intelligence application for your daily Phage Research

**Description:**

https://phage.ai/accounts/login/?next=/
https://github.com/phageaisa/phageai

PhageAI is an application that simultaneously represents a repository of knowledge of bacteriophages and a tool to analyse genomes with Artificial Intelligence support.

**Framework modules**

Set of methods related with:

* `lifecycle` - bacteriophage lifecycle prediction:

  * `.predict(fasta_path)` - return bacteriophage lifecycle prediction class (Virulent, Temperate or Chronic) with probability (%); 
* `taxonomy` - bacteriophage taxonomy order, family and genus prediction (TBA);  
* `topology` - bacteriophage genome topology prediction (TBA);  
* `repository` - set of methods related with PhageAI bacteriophage repository:

  * `.get_record(value)` - return dict with Bacteriophage meta-data
  * `.get_top10_similar_phages(value)` - return list of dicts contained top-10 most similar bacteriophages

------

Machine Learning algorithms can process enormous amounts of data in relatively short time in order to find connections and dependencies that are unobvious for human beings. Correctly designed applications based on AI are able to vastly improve and speed up the work of the domain experts.

Models based on DNA contextual vectorization and Deep Neural Networks are particularly effective when it comes to analysis of genomic data. The system that we propose aims to use the phages sequences uploaded to the database to build a model which is able to predict if a bacteriophage is virulent, temperate or chronic with a high probability.

One of the key system modules is the bacteriophages repository with a clean web interface that allows to browse, upload and share data with other users. The gathered knowledge about the bacteriophages is not only valuable on its own but also because of the ability to train the ever-improving Machine Learning models.

Detection of virulent or temperate features is only one of the first tasks that can be solved with Artificial Intelligence. The combination of Biology, Natural Language Processing and Machine Learning allows us to create algorithms for genomic data processing that could eventually turn out to be effective in a wide range of problems with focus on classification and information extraced from DNA.


### PHACTS

**Description:**

https://pubmed.ncbi.nlm.nih.gov/22238260/
https://edwards.sdsu.edu/PHACTS/
https://edwards.sdsu.edu/PHACTS/PHACTS-0.3.tar.gz

**Abstract**
*Motivation*: Bacteriophages have two distinct lifestyles: virulent and temperate. The virulent lifestyle has many implications for phage therapy, genomics and microbiology. Determining which lifestyle a newly sequenced phage falls into is currently determined using standard culturing techniques. Such laboratory work is not only costly and time consuming, but also cannot be used on phage genomes constructed from environmental sequencing. Therefore, a computational method that utilizes the sequence data of phage genomes is needed.

*Results*: Phage Classification Tool Set (PHACTS) utilizes a novel similarity algorithm and a supervised Random Forest classifier to make a prediction whether the lifestyle of a phage, described by its proteome, is virulent or temperate. The similarity algorithm creates a training set from phages with known lifestyles and along with the lifestyle annotation, trains a Random Forest to classify the lifestyle of a phage. PHACTS predictions are shown to have a 99% precision rate.

*Availability and implementation*: PHACTS was implemented in the PERL programming language and utilizes the FASTA program (Pearson and Lipman, 1988) and the R programming language library 'Random Forest' (Liaw and Weiner, 2010). The PHACTS software is open source and is available as downloadable stand-alone version or can be accessed online as a user-friendly web interface. The source code, help files and online version are available at http://www.phantome.org/PHACTS/.



### phap - Phage Host Analysis Pipeline

**Description:**

https://github.com/MGXlab/phap
PHAP wraps the execution of various phage-host prediction tools.

Overview
Features
Uses Singularity containers for the execution of all tools.

When possible (i.e. the image is not larger than a few Gs), tools and their dependencies are bundled in the same container. This means you do not need to get models or any other external databases, unless otherwise specified.

Intermediate processing steps are handled by Conda environments, to ensure smooth and reproducible execution.

Outputs the Last Common Ancestor of all tools, per contig, based on the predicted taxonomy.

* HTP https://github.com/wojciech-galan/viruses_classifier
* RaFAh https://sourceforge.net/projects/rafah/
* vHuLK https://github.com/LaboratorioBioinformatica/vHULK
* VirHostMatcher-Net https://github.com/WeiliWw/VirHostMatcher-Net
* WIsH https://github.com/soedinglab/WIsH



### PHERI - Phage Host Exploration pipeline

**Description:**

https://www.biorxiv.org/content/10.1101/2020.05.13.093773v3.full

The solution to this problem may be to use a bioinformatic approach in the form of prediction software capable of determining a bacterial host based on the phage whole-genome sequence. The result of our research is the machine learning algorithm based tool called PHERI. PHERI predicts suitable bacterial host genus for purification of individual viruses from different samples. Besides, it can identify and highlight protein sequences that are important for host selection. PHERI is available at 
https://hub.docker.com/repository/docker/andynet/pheri. 

The source code for the model training is available at
https://github.com/andynet/pheri_preprocessing

, and the source code for the tool is available at
https://github.com/andynet/pheri.


### Phirbo - A tool to predict prokaryotic hosts for phage (meta)genomic sequences

**Description:**

https://github.com/aziele/phirbo

Phirbo links phage to host sequences through other intermediate sequences that are potentially homologous to both phage and host sequences.

To link phage (P) to host (H) sequence through intermediate sequences, phage and host sequences need to be used as queries in two separate sequence similarity searches (e.g., BLAST) against the same reference database of prokaryotic genomes (D). One BLAST search is performed for phage query (P) and the other for host query (H). The two lists of BLAST results, P → D and H → D, contain prokaryotic genomes ordered by decreasing score. To avoid a taxonomic bias due to multiple genomes of the same prokaryote species (e.g., Escherichia coli), prokaryotic species can be ranked according to their first appearance in the BLAST list. In this way, both ranked lists represent phage and host profiles consisting of the ranks of top-score prokaryotic species.

Phirbo estimates the phage-host relationship by comparing the content and order between phage and host ranked lists using Rank-Biased Overlap (RBO) measure. Briefly, RBO fosters comparison of ranked lists of different lengths with heavier weights for matching the higher-ranking items. RBO ranges between 0 and 1, where 0 means that the lists are disjoint (have no items in common) and 1 means that the lists are identical in content and order.


### PHIST - Phage-Host Interaction Search Tool

**Description:**

https://github.com/refresh-bio/PHIST
https://www.biorxiv.org/content/10.1101/2021.09.06.459169v1

A tool to predict prokaryotic hosts for phage (meta)genomic sequences. PHIST links viruses to hosts based on the number of k-mers shared between their sequences.


### rafah - Random Forest Assignment of Hosts

**Description:**

https://www.sciencedirect.com/science/article/pii/S2666389921001008
https://sourceforge.net/projects/rafah/

One fundamental question when trying to describe viruses of Bacteria and Archaea is: Which host do they infect? To tackle this issue we developed a machine-learning approach named Random Forest Assignment of Hosts (RaFAH), which outperformed other methods for virus-host prediction. Our rationale was that the machine could learn the associations between genes and hosts much more efficiently than a human, while also using the information contained in the hypothetical proteins. Random forest models were built using the Ranger⁠ package in R⁠.


### VirHostMatcher-Net

**Description:**

https://github.com/WeiliWw/VirHostMatcher-Net
https://academic.oup.com/nargab/article/2/2/lqaa044/5861484

Metagenomic sequencing has greatly enhanced the discovery of viral genomic sequences; however, it remains challenging to identify the host(s) of these new viruses. We developed VirHostMatcher-Net, a flexible, network-based, Markov random field framework for predicting virus–prokaryote interactions using multiple, integrated features: CRISPR sequences and alignment-free similarity measures (⁠s∗2 and WIsH). Evaluation of this method on a benchmark set of 1462 known virus–prokaryote pairs yielded host prediction accuracy of 59% and 86% at the genus and phylum levels, representing 16–27% and 6–10% improvement, respectively, over previous single-feature prediction approaches. We applied our host prediction tool to crAssphage, a human gut phage, and two metagenomic virus datasets: marine viruses and viral contigs recovered from globally distributed, diverse habitats. Host predictions were frequently consistent with those of previous studies, but more importantly, this new tool made many more confident predictions than previous tools, up to nearly 3-fold more (n > 27 000), greatly expanding the diversity of known virus–host interactions.


### vHULK

**Description:**

https://github.com/LaboratorioBioinformatica/vHULK
https://www.biorxiv.org/content/10.1101/2020.12.06.413476v1.full

**Phage Host Prediction using high level features and neural networks**

Metagenomics and sequencing techniques have greatly improved in these last five years and, as a consequence, the amount of data from microbial communities is astronomic. An import part of the microbial community are phages, which have their own ecological roles in the environment. Besides that, they have also been given a possible human relevant (clinical) role as terminators of multidrug resistant bacterial infections. A lot of basic research still need to be done in the Phage therapy field, and part of this research involves gathering knowledge from new phages present in the environment as well as about their relationship with clinical relevant bacterial pathogens.

Having this scenario in mind, we have developed vHULK. A user-friendly tool for prediction of phage hosts given their complete or partial genome in FASTA format. Our tool outputs an ensemble prediction at the genus or species level based on scores of four different neural network models. Each model was trained with more than 4,000 genomes whose phage-host relationship was known. v.HULK also outputs a mesure of entropy for each final prediction, which we have demonstrated to be correlated with prediction's accuracy. The user might understand this value as additional information of how certain v.HULK is about a particular prediction. We also suspect that phages with higher entropy values may have a broad host-range. But that hypothesis is to be tested later. Accuracy results in test datasets were >99% for predictions at the genus level and >98% at the species level. vHULK currently supports predictions for 52 different prokaryotic host species and 61 different genera.


### WIsH - who is the host? Predicting prokaryotic hosts from metagenomic phage contigs

**Description:**

https://doi.org/10.1093/bioinformatics/btx383

https://github.com/soedinglab/WIsH

WIsH can identify bacterial hosts from metagenomic data, keeping good accuracy even on 
smaller contigs.

WIsH predicts prokaryotic hosts of phages from their genomic sequences. It achieves 63% mean accuracy when predicting the host genus among 20 genera for 3 kbp-long phage contigs. Over the best current tool, WisH shows much improved accuracy on phage sequences of a few kbp length and runs hundreds of times faster, making it suited for metagenomics studies.

### Digital Phagogram