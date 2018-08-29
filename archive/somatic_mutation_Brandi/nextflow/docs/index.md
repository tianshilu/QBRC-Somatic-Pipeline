# Astrocyte Germline Variant Calling Workflow Package

This workflow carries out a Germline Exome Analysis pipeline, including the integration of variants from various callers and basic annotation.

    1) Reads are trimmed by quality (threshold of 25) and adapter sequence at the ends and reads shorter than 35bps are removed.

    2) Trimmed Fastq files are aligned to the selected reference genome using BWA-MEM (Li et al 2013)

    3) Subsequent processing will be performed by Samtools and Picard to ensure proper file formatting and mark duplicates.
 
    4) Alignments are then recalibrated and realigned using GATK (DePristo et al 2011;McKenna et al 2010)

    5) To detect genome germline variants, GATK2 (DePristo et al 2011, McKenna et al 2010), Platypus (Rimmer et al 2014), Samtools version 1.3 and FreeBayes version 0.9.7 (Garrison and Marth 2012) are used. 

    6) Integration of predicted SNPs and INDELs from these algorithms is performed using BAYSIC (Cantarel et al 2014).

    7) Effect of SNPs and INDELs on genes is predicted using snpEff (Cingolani et al 2012) using the GRCh38.82 database. Allele frequency in the general population is determined by comparison to ExAC (The ExAC Consortium 2015). Additionally, discovered variants are annotated using SnpSift (Cingolani et al 2012) the dbSNP, COSMIC (Forbes et al. 2009), CLINVAR (Landrum et al 2014), GWAS Catalog (Welter et al 2014) and dbNSFP (Liu et al 2011) databases.

##Workflow Parameters

    fastq - Choose one or more RNASeq read files to process.

    genome - Choose a genomic reference (genome).

    pairs - Choose if pair-ended or single-end sequences

    capturebed - Choose the bedfile that corresponds to the capture kit used for the sequencing.  This is used to calculate average sequence coverage and filter SNPs

    design - This file matches the fastq files to data about the sample

 The following columns are necessary, must be named as in template and can be in any order:

    SampleID
        This ID should match the name in the fastq file ie S0001.R1.fastq.gz the sample ID is S0001
    SampleName
        This ID can be the identifier of the researcher or clinician
    SubjectID
        Used in order to link samples from the same patient
    Phenotype
	2= Case or Diseaes Phenotype, 1= Healthy Control
    Gender
	1=male, 2=female
    FullPathToFqR1
	Name of the fastq file R1
    FullPathToFqR2
	Name of the fastq file R2

There are some optional columns that might help with the analysis:
      SequenceRun
      Organism
      FamilyID
      CellPopulation
      Treatment
      GeneticFeature (WT or KO)
      Race
      Ethnicity
      Age
      

### Test Data


### Credits
This example worklow is derived from original scripts kindly contributed by the Bioinformatic Core Facility (BICF), Department of Bioinformatics and Clinical Sciences.

### References

    Andy Rimmer, Hang Phan, Iain Mathieson, Zamin Iqbal, Stephen R. F. Twigg, WGS500 Consortium, Andrew O. M. Wilkie, Gil McVean, Gerton Lunter. Integrating mapping-, assembly- and haplotype-based approaches for calling variants in clinical sequencing applications. Nature Genetics (2014) doi:10.1038/ng.3036
    Bernstein, B. E., Birney, E., Dunham, I., Green, E. D., Gunter, C., & Snyder, M. (2012). An integrated encyclopedia of DNA elements in the human genome. Nature, 489, 57–74. doi:10.1038/nature11247
    Cantarel, B. L., Weaver, D., McNeill, N., Zhang, J., Mackey, A. J., & Reese, J. (2014). BAYSIC: a Bayesian method for combining sets of genome variants with improved specificity and sensitivity. BMC Bioinformatics, 15, 104. doi:10.1186/1471-2105-15-104
    Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L., ? Ruden, D. M. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly. doi:10.4161/fly.19695
    Cingolani, P., Patel, V. M., Coon, M., Nguyen, T., Land, S. J., Ruden, D. M., & Lu, X. (2012). Using Drosophila melanogaster as a Model for Genotoxic Chemical Mutational Studies with a New Program, SnpSift. Frontiers in Genetics. doi:10.3389/fgene.2012.00035
    Challis, D., Yu, J., Evani, U. S., Jackson, A. R., Paithankar, S., Coarfa, C., ? Yu, F. (2012). An integrative variant analysis suite for whole exome next-generation sequencing data. BMC Bioinformatics. doi:10.1186/1471-2105-13-8
    DePristo, M. A., Banks, E., Poplin, R., Garimella, K. V, Maguire, J. R., Hartl, C., ? Daly, M. J. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics, 43, 491–498. doi:10.1038/ng.806
    Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) [01 (01, 2016) accessed].
    Forbes, S. A., Tang, G., Bindal, N., Bamford, S., Dawson, E., Cole, C., ? Futreal, P. A. (2009). COSMIC (the Catalogue of Somatic Mutations In Cancer): A resource to investigate acquired mutations in human cancer. Nucleic Acids Research, 38. doi:10.1093/nar/gkp995
    Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012
    Hansen NF, Gartner JJ, Mei L, Samuels Y, Mullikin JC. Shimmer: detection of genetic alterations in tumors using next-generation sequence data. Bioinformatics. 2013 Jun 15;29(12):1498-503. doi: 10.1093/bioinformatics/btt183. Epub 2013 Apr 24. PubMed PMID: 23620360; PubMed Central PMCID: PMC3673219.
    Kim S, Jeong K, Bhutani K, Lee J, Patel A, Scott E, Nam H, Lee H, Gleeson JG, Bafna V. Virmid: accurate detection of somatic mutations with sample impurity inference. Genome Biol. 2013 Aug 29;14(8):R90. doi: 10.1186/gb-2013-14-8-r90. PubMed PMID: 23987214; PubMed Central PMCID: PMC4054681.
    Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Res. 2012 Mar;22(3):568-76. doi: 10.1101/gr.129684.111. Epub 2012 Feb 2. PubMed PMID: 22300766; PubMed Central PMCID: PMC3290792.
    Landrum, M. J., Lee, J. M., Riley, G. R., Jang, W., Rubinstein, W. S., Church, D. M., & Maglott, D. R. (2014). ClinVar: Public archive of relationships among sequence variation and human phenotype. Nucleic Acids Research, 42. doi:10.1093/nar/gkt1113
    Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv Preprint arXiv, 00, 3. doi:arXiv:1303.3997 [q-bio.GN]
    Liu X, Jian X, Boerwinkle E. dbNSFP: a lightweight database of human nonsynonymous SNPs and their functional predictions. Hum Mutat. 2011 Aug;32(8):894-9. doi: 10.1002/humu.21517. PubMed PMID: 21520341
    McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ? DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20, 1297–1303. doi:10.1101/gr.107524.110
    The 1000 Genome Consortium. An integrated map of genetic variation from 1,092 human genomes. (2012). Nature, 491(7422), 56–65. Retrieved from http://dx.doi.org/10.1038/nature11632.
    Saunders CT, Wong WS, Swamy S, Becq J, Murray LJ, Cheetham RK. Strelka: accurate somatic small-variant calling from sequenced tumor-normal sample pairs. Bioinformatics. 2012 Jul 15;28(14):1811-7. doi: 10.1093/bioinformatics/bts271. Epub 2012 May 10. PubMed PMID: 22581179.
    Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H. The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006.
