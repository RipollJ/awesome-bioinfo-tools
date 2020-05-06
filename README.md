# Awesome Bioinfo-tools

A curated list of awesome Bioinformatics tools.

----------------------------------------------------------------------------
Table of contents
- [Awesome Bioinfo-tools](#awesome-bioinfo-tools)
  * [Awesome existing topics related to bioinformatics](#awesome-existing-topics-related-to-bioinformatics)
  * [Suite tools](#suite-tools)
  * [Quality analysis & trimming tools](#quality-analysis---trimming-tools)
    + [quality analysis checking](#quality-analysis-checking)
    + [trimming](#trimming)
    + [read merger](#read-merger)
    + [demultiplexing](#demultiplexing)
  * [Multiviewer](#multiviewer)
  * [Mapping tools](#mapping-tools)
    + [mapping](#mapping)
    + [splice-aligner](#splice-aligner)
  * [Assembly tools](#assembly-tools)
    + [Genome & Transcriptome assembly](#genome---transcriptome-assembly)
    + [Metagenome & Metatranscriptome assembly](#metagenome---metatranscriptome-assembly)
    + [Viewers](#viewers)
    + [Correction tools](#correction-tools)
  * [Variant calling & alternative splicing tools](#variant-calling---alternative-splicing-tools)
    + [variant calling](#variant-calling)
    + [Motif discovery](#motif-discovery)
    + [Peak calling](#peak-calling)
    + [Learning tools](#learning-tools)
    + [Correction tools](#correction-tools-1)
  * [Counting tools](#counting-tools)
  * [Differential analysis tools](#differential-analysis-tools)
    + [RNA-seq](#rna-seq)
    + [Alternative-splicing](#alternative-splicing)
    + [RIBO-seq](#ribo-seq)
  * [Others](#others)
    + [Exploration tools for RNA-seq and RIBO-seq](#exploration-tools-for-rna-seq-and-ribo-seq)
    + [Network visualisation](#network-visualisation)
    + [Interaction visualisation](#interaction-visualisation)
    + [Annotations & pathway databases](#annotations---pathway-databases)
  * [Bioinformatic analysis informations (Wikipedia links)](#bioinformatic-analysis-informations--wikipedia-links-)
    + [Metagenomic](#metagenomic)
    + [Metatranscriptomic](#metatranscriptomic)
    + [Metabarcoding](#metabarcoding)
    + [Alternative-splicing](#alternative-splicing-1)
    + [Ribo-seq](#ribo-seq)
    + [Merip-seq](#merip-seq)
    + [mi-CLIP](#mi-clip)
    + [Proteomics](#proteomics)
    + [MASS-SPEC](#mass-spec)
  * [Specific workflow or platform](#specific-workflow-or-platform)

----------------------------------------------------------------------------
## Awesome existing topics related to bioinformatics
- [awesome-bioinformatics](https://github.com/danielecook/Awesome-Bioinformatics): some informations on Bioinformatics.
- [awesome-alternative-splicing](https://github.com/HussainAther/awesome-alternative-splicing): some programs for alternative splicing analysis.
- [awesome-pipeline](https://github.com/pditommaso/awesome-pipeline): found the best pipeline workflow for your applications.
- [awesome-deep-learning](https://github.com/ChristosChristofidis/awesome-deep-learning): deep learning informations.

[<small>[top↑]</small>](#)
----------------------------------------------------------------------------
## Suite tools
- [BBtools](https://jgi.doe.gov/data-and-tools/bbtools/): BBTools is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA sequence data. 
- [samtools](https://github.com/samtools/samtools): The original samtools package has been split into three separate but tightly coordinated projects:
    - htslib: C-library for handling high-throughput sequencing data
    - samtools: mpileup and other tools for handling SAM, BAM, CRAM
    - bcftools: calling and other tools for handling VCF, BCF
- [GATK](https://gatk.broadinstitute.org/hc/en-us): A genomic analysis toolkit focused on variant discovery.
- [EA-Utils](https://expressionanalysis.github.io/ea-utils/): Command-line tools for processing biological sequencing data. Barcode demultiplexing, adapter trimming, etc. Primarily written to support an Illumina based pipeline - but should work with any FASTQs.

[<small>[top↑]</small>](#)
----------------------------------------------------------------------------
## Quality analysis & trimming tools
### quality analysis checking
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): A quality control tool for high throughput sequence data.
- [FastQ Screen](http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html): FastQ Screen is a simple application which allows you to search a large sequence dataset against a panel of different genomes to determine from where the sequences in your data originate.

### trimming
- [Sickle](https://github.com/najoshi/sickle): A windowed adaptive trimming tool for FASTQ files using quality.
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html): Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
- [bbduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/): “Duk” stands for Decontamination Using Kmers. BBDuk was developed to combine most common data-quality-related trimming, filtering, and masking operations into a single high-performance tool. 
- [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/): A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries.
- [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic): A flexible read trimming tool for Illumina NGS data.
- [Sortmerna](https://bioinfo.lifl.fr/RNA/sortmerna/): Fast filtering, mapping and OTU picking.

### read merger
- [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html): PEAR is an ultrafast, memory-efficient and highly accurate pair-end read merger. 
- [Fastq-join](https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqJoin.md): Joins two paired-end reads on the overlapping ends.
- [Seq-prep](https://github.com/jstjohn/SeqPrep): SeqPrep is a program to merge paired end Illumina reads that are overlapping into a single longer read.
- [FLASH](http://ccb.jhu.edu/software/FLASH/): Fast Length Adjustment of SHort reads

### demultiplexing
- [fastq-multx](https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMultx.md): The goal of this program is to make it easier to demultiplex possibly paired-end sequences, and also to allow the "guessing" of barcode sets based on master lists of barcoding protocols (fluidigm, truseq, etc.)
- [UMI-tools](https://github.com/CGATOxford/UMI-tools): This repository contains tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs) and single cell RNA-Seq cell barcodes.
- [NextGenSeqUtils](https://nextjournal.com/Murrell-Lab/demultiplexing-with-custom-barcoded-primers): Notebook for demultiplexing with custom barcoded primers.

[<small>[top↑]</small>](#)
----------------------------------------------------------------------------
## Multiviewer
- [MultiQC](https://multiqc.info/): Aggregate results from bioinformatics analyses across many samples into a single report.

[<small>[top↑]</small>](#)
----------------------------------------------------------------------------
## Mapping tools
### mapping
- [BWA](https://github.com/lh3/bwa): BWA is a software package for mapping DNA sequences against a large reference genome, such as the human genome. 
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.
- [PANDASEQ](https://github.com/neufeld/pandaseq): PANDASEQ is a program to align Illumina reads, optionally with PCR primers embedded in the sequence, and reconstruct an overlapping sequence.

### splice-aligner
- [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml): TopHat is a fast splice junction mapper for RNA-Seq reads.
- [STAR](https://github.com/alexdobin/STAR): Spliced Transcripts Alignment to a Reference.
- [CRAC](http://crac.gforge.inria.fr/): RNA-Seq mapping software that include the discovery of transcriptomic and genomic variants like splice junction, chimeric junction, SNVs, Indels in a single analysis step using a built-in error detection method enabling high precison and sensitivity.

[<small>[top↑]</small>](#)
---------------------------------------------------------------------------
## Assembly tools
### Genome & Transcriptome assembly
- [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/): Sequence assembler for very short reads
- [SPAdes](http://cab.spbu.ru/files/release3.12.0/manual.html): SPAdes – St. Petersburg genome assembler – is an assembly toolkit containing various assembly pipelines.
- [Minia](https://github.com/GATB/minia): Minia is a short-read assembler based on a de Bruijn graph, capable of assembling a human genome on a desktop computer in a day.

### Metagenome & Metatranscriptome assembly
- [Megahit](https://github.com/voutcn/megahit): An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph
- [MetaVelvetSL](http://metavelvet.dna.bio.keio.ac.jp/MSL.html): An extension of Velvet assembler to de novo metagenomic assembly
- [MetaSPADES](https://kbase.us/applist/apps/kb_SPAdes/run_metaSPAdes/release?gclid=CjwKCAjwwMn1BRAUEiwAZ_jnEkxHhKIJi9wEvB5wcMs818fq2fFTmdO_-2Dz3EP015_QPbgdBO73tRoCdxoQAvD_BwE): Assemble metagenomic reads using the SPAdes assembler.
- [Minia for metagenome](https://github.com/GATB/gatb-minia-pipeline): GATB-Minia-Pipeline is a de novo assembly pipeline for Illumina data. It can assemble genomes and metagenomes.

### Viewers
- [Bandage](https://rrwick.github.io/Bandage/): Bandage is a program for visualising de novo assembly graphs.
- [IGV](http://software.broadinstitute.org/software/igv/): visualization tool for interactive exploration of large, integrated genomic datasets.

### Correction tools
- [rnaQUAST](http://cab.spbu.ru/software/rnaquast/): rnaQUAST is a software designed for quality evaluation and assessment of de novo transcriptome assemblies. 

[<small>[top↑]</small>](#)
---------------------------------------------------------------------------
## Variant calling & alternative splicing tools
### variant calling
- [VarScan](http://varscan.sourceforge.net/): variant detection in massively parallel sequencing data.
- [KisSplice](http://kissplice.prabi.fr/): A local transcriptome assembler for SNPs, indels and AS events
- [Farline](http://fasterdb.ens-lyon.fr/tools.php): FaRLine is a pipeline to analyse the alternative splicing.
- [SplAdder](https://github.com/ratschlab/spladder): SplAdder, short for Splicing Adder, a toolbox for alternative splicing analysis based on RNA-Seq alignment data.
- [Whippet](https://github.com/timbitz/Whippet.jl): Efficient and Accurate Quantitative Profiling of Alternative Splicing Patterns of Any Complexity on a Laptop.
- [freebayes](https://github.com/ekg/freebayes/blob/master/README.md): freebayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

### Motif discovery
- [MaxENTScan](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html): MaxEntScan is based on the approach for modeling the sequences of short sequence motifs such as those involved in RNA splicing which simultaneously accounts for non-adjacent as well as adjacent dependencies between positions.

### Peak calling
- [MACS2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html): computational method used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing experiment.
- [m6a Viewer](http://dna2.leeds.ac.uk/m6a/): m6a Viewer is a cross-platform java application for detecting and visualising peaks in ME-RIP/ m6a-seq data. 

### Learning tools
- [DeepVariant](https://github.com/google/deepvariant): DeepVariant is an analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data. 
- [LaBrachoR](https://github.com/jpaggi/labranchor): LaBranchoR uses a LSTM network built with keras to predict the position of RNA splicing branchpoints relative to a three prime splice site.
- [SpliceAI](https://github.com/Illumina/SpliceAI): A deep learning-based tool to identify splice variants
- [SpliceAI-wrapper](https://github.com/bihealth/spliceai-wrapper): SpliceAI Wrapper, is an attempt to use caching for reducing the number of required predictions. Please note that the authors of SpliceAI Wrapper are unrelated to the authors of SpliceAI.

### Correction tools
- [Portcullis](https://github.com/TGAC/portcullis): Portcullis stands for PORTable CULLing of Invalid Splice junctions from pre-aligned RNA-seq data. 

[<small>[top↑]</small>](#)
-------------------------------------------------------------------------
## Counting tools
- [FeatureCounts](http://bioinf.wehi.edu.au/featureCounts/): counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations.
- [Kallisto](https://pachterlab.github.io/kallisto/about): kallisto is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads. 
- [HTSeqCount](https://htseq.readthedocs.io/en/master/): Analysing high-throughput sequencing data with Python
- [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#top): Transcript assembly and quantification for RNA-Seq

[<small>[top↑]</small>](#)
-------------------------------------------------------------------------
## Differential analysis tools
### RNA-seq
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html): Differential gene expression analysis based on the negative binomial distribution.
- [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html): Empirical Analysis of Digital Gene Expression Data in R.
- [NBAMSeq](https://bioconductor.org/packages/release/bioc/html/NBAMSeq.html): NBAMSeq is a Bioconductor package for differential expression analysis based on negative binomial additive model.
- [NOISeq](http://bioconductor.org/packages/release/bioc/html/NOISeq.html): Exploratory analysis and differential expression for RNA-seq data
- [Sleuth](https://pachterlab.github.io/sleuth/about): sleuth is a program for analysis of RNA-Seq experiments for which transcript abundances have been quantified with kallisto. 

### Alternative-splicing
- [DEX-seq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html): Inference of differential exon usage in RNA-Seq.
- [KissDE](https://bioconductor.org/packages/release/bioc/html/kissDE.html): Retrieves Condition-Specific Variants in RNA-Seq Data.

### RIBO-seq
- [Xtail](https://github.com/xryanglab/xtail): Genome-wide assessment of differential translations with ribosome profiling data.
- [Anota2Seq](https://bioconductor.org/packages/release/bioc/html/anota2seq.html): Generally applicable transcriptome-wide analysis of translational efficiency using anota2seq.

[<small>[top↑]</small>](#)
------------------------------------------------------------------------
## Others
### Exploration tools for RNA-seq and RIBO-seq
- [RNA-Ribo Explorer (RRE)](http://www.atgc-montpellier.fr/rre/): RRE is an interactive, stand-alone, and graphical software for analysing, viewing and mining both transcriptome (typically RNA-seq) and translatome (typically Ribosome profiling or Ribo-seq) datasets.
- [IGET](https://iget.c2b2.columbia.edu/): The Integrated Genomics Exploration Tools (IGET) website provides online access to a suite of tools for exploring biological pathways and DNA/RNA/protein regulatory elements associated with large-scale gene expression and protein behavior dynamics.

### Network visualisation
- [Gephi](https://gephi.org/): visualization and exploration software for all kinds of graphs and networks. 
- [Cytoscape](https://cytoscape.org/): visualization of complex networks and integrating these with any type of attribute data. 

### Interaction visualisation
- [String](https://string-db.org/): Protein-Protein Interaction Networks Functional Enrichment Analysis

### Annotations & pathway databases
- [Gene Ontology](http://geneontology.org/): The Gene Ontology (GO) knowledgebase is the world’s largest source of information on the functions of genes.
- [gProfiler](https://biit.cs.ut.ee/gprofiler/gost): g:Profiler is a public web server for characterising and manipulating gene lists. 
- [KEGG](https://www.genome.jp/kegg/): KEGG is a database resource for understanding high-level functions and utilities of the biological system, such as the cell, the organism and the ecosystem, from molecular-level information, especially large-scale molecular datasets generated by genome sequencing and other high-throughput experimental technologies.
- [DAVID](https://david.ncifcrf.gov/): Database for Annotation, Visualization and Integrated Discovery (DAVID).
- [PANTHER](http://www.pantherdb.org/): Protein ANalysis THrough Evolutionary Relationships.
- [RNAcentral](https://rnacentral.org/): The non-coding RNA sequence database

[<small>[top↑]</small>](#)
-------------------------------------------------------------------------
## Bioinformatic analysis informations (Wikipedia links)
### Metagenomic
- [Metagenomics](https://en.wikipedia.org/wiki/Metagenomics)

### Metatranscriptomic
- [Metatranscriptomics](https://en.wikipedia.org/wiki/Metatranscriptomics)

### Metabarcoding
- [Barcoding](https://fr.wikipedia.org/wiki/Barcoding_mol%C3%A9culaire)

### Alternative-splicing
- [Alternative-splicing](https://en.wikipedia.org/wiki/Alternative_splicing)

### Ribo-seq
- [Ribosome profiling](https://en.wikipedia.org/wiki/Ribosome_profiling)

### Merip-seq
- [Merip-seq](https://en.wikipedia.org/wiki/MeRIPseq)

### mi-CLIP
- [Illumina mi-CLIP-m6a](https://emea.illumina.com/science/sequencing-method-explorer/kits-and-arrays/miclip-m6a.html)

### Proteomics
- [Proteomics](https://en.wikipedia.org/wiki/Proteomics)

### MASS-SPEC
- [Mass spectrometry in general](https://en.wikipedia.org/wiki/Mass_spectrometry)
- [Mass-spec for DNA](https://www.ncbi.nlm.nih.gov/pubmed/15829234)
- [Mass-spec for RNA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2442014/)

[<small>[top↑]</small>](#)
----------------------------------------------------------------------------
## Specific workflow or platform
- [KisSplice](http://kissplice.prabi.fr/training/): Training alternative splicing analysis with KisSplice & suite tools.
- [QIIME2](https://qiime2.org/): QIIME 2™ is a next-generation microbiome bioinformatics platform that is extensible, free, open source, and community developed.
- [Mothur](https://mothur.org/): This project seeks to develop a single piece of open-source, expandable software to fill the bioinformatics needs of the microbial ecology community.
- [Vsearch](https://github.com/torognes/vsearch): Open source tool for metagenomics.

[<small>[top↑]</small>](#)