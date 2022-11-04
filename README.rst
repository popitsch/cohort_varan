   ::

    .--------------------------.
    | c o h o r t _ v a r a n  |\\
    |   __{    }__ __{    }__  | |
    |  || ( .. ) ||  ( .. ) || | |
    |  __  \__/  __   \__/  __ | |
    | /..\  ''  /..\   ''  /..\| |
    .--------------------------. |
     \__________v_1.0___________\|
 
Introduction
============
 
Cohort_varan consists of two python pipelines for the comprehensive annotation of single VCF files (annotate_vcf)
and sets of VCF files (annotate_cohort).

Overall usage scenario for cohort analyses:

* Annotate small (SNV and INDEL) and large (SV) VCF files with annotate_vcf
* Extract data per pedigree and integrate with various datasources with annotate_cohort
* Preprocess data with `Variant Explorer R script <https://github.com/edg1983/Variant_explorer/blob/master/preprocessing/Prepare_Rdata_object.R>`_
* Load and analyse data with `Variant Explorer <https://github.com/edg1983/Variant_explorer>`_
 

Installation
============

* install python 3.6+ and clone this repository
* install the following 3rd party tools

    - `htslib <http://www.htslib.org/>`_
    - `SnpEff and SnpSift <http://pcingola.github.io/SnpEff/>`_
    - `vcfanno <https://github.com/brentp/vcfanno>`_
    - `vcftools <https://vcftools.github.io/index.html>`_ (for vcf-sort support; only used for platypus VCF files)
 
annotate_vcf
============
 
This is a python pipeline for the comprehensive annotation of VCF files using SnpEff/SnpSift and vcfanno.
Briefly, the pipeline consists of the following steps:

* Links to VCF files are read from a sample sheet and iteratively processed. 
* VCF files are pre-filtered to

    - remove non-PASS variants (unless the 'includeNonPass' option is set)
    - contain only variants for the passed chromosomes
* VCF files are preprocessed in a variant-caller specific manner:

    - platypus VCF files: files are sorted with vcf-sort
    - deepvariant/gatk VCF files: entries in the ID column are moved to an info field
* VCF files are annotated with snpEFF
* LUA files for the efficient annotation with vcfanno are created for the configured resources 
* VCF files are annotated with vcfanno
* VCF files are post-filtered with snpSift using a configurable filter string
* Configured fields from the final VCF (including sample-specific values, e.g., genotypes) are extracted and written to a TSV file.
   

Example call
------------

  ::

     python cohort_varan/cohort_varan.py annotate_vcf \
        -i samples.tsv \
        -o annotated_vcf_files \
        -c config.json \
        --snpEffConfig <path>/snpEff/snpEff.config \
        --threads 2


Sample sheet 
------------
Sample sheets are TSV files with 4 columns:

* id: used as output filename prefix
* genome: string that must match the respective entry in the VCF config file, e.g., "GRCh38" (see below)
* source: 'platypus'|'deepvariant'|'gatk'. Used for caller-specific VCF preprocessing.
* vcf_file: path to VCF file

Example:

  ::

    #PID\tgenome\tsource\tvcf
    test\tGRCh38\tdeepvariant\ttest.vcf

This configures a single VCF file 'test.vcf' that will be annotated with the config entries from the 'GRCh38' block. The
ID entries in the input file will be moved to a field 'old_deepvariant_ids'.

Config file
-----------
JSON file with the following structure (cf. example below):

* ``prefix``: in this section, variables can be defined that will be replaced throughout the remaining config file. 

    - Example: "prefix" : { "@REF": "/path_prefix" } 
* ``ref``: reference genome definitions. The contained section names must match the entries in the sample sheet. subsections: 

    - FASTA (reference genomme fasta file)
    - snpEff (snpEff annotation db version)
    - chr (list of considered chromosome names
    - filter (snpEff filter string)
* ``roi``: regions of interest annotation BED files. Must have subsections for each configured reference genome.
* ``af``: population allele frequency VCF files.  Must have subsections for each configured reference genome and beyond this 2 subsections "global" (global AF annotations) and "pop" (population specific AF annotations). Each entry has 2 parts: VCF file and comma-separated list of info fields to be loaded from these vcf files
* ``anno``: general annotation files. Must have subsections for each configured reference genome. Each entry has multiple parts, depending on file type:

   - BED files: [file_path]
   - VCF files: [file_path, comma-separated_list_of_info_field_names,comma-separated_list_of_summarization_methods see vcfanno docs; if omitted, 'self' will be used)]                                
   - TSV files: [file_path,comma-separated_list_of_column_indices,comma-separated_list_of_variable_names (if omitted, the TSV column names will be used), comma-separated_list_of_summarization_methods (see vcfanno docs; if omitted, 'self' will be used)]
* ``known``: VCF files with 'known' variants.  Must have subsections for each configured reference genome. The respective IDs will be written to the annotated VCF ID section. It is possible to configure a prefix string for these IDs
* ``output``: List of output fields for the created TSV file. There must be a 'fields' section containing list entries with a column name and an optional operator: 

    - 'max': split by comma, replace '.' values with nan (will be ignored) and select maximum value
    - 'min': split by comma, replace '.' values with nan (will be ignored) and select minimum value
    - none: use value as is.
* ``tools``: Optional section for providing custom paths for the following 3rd party tools:

    - vcf-sort
    - snpSift
    - snpEff
    - vcfanno

    If omitted, the pipeline will try to call the tool directly by name.
* ``linux_temp_dir``: optional, for configuring an alternative TEMP dir. 
    

Example JSON config file
------------------------

`<annotate_vcf.example_config.json>`_


This file configures a single reference genome (GRCh38) and uses snpEff database 'GRCh38.99'. Only 2 chromosomes and unfiltered variants with a quality >= 10 will be included.Regions of interest will be read from a BED file and the respective path will be prefixed by the configured "@REF" prefix. Global and population-specific allele frequencies will be read from the configured info fields ('AF' and 'AF_EUR') in the configured VCF files. CADD scores (raw and phred-scaled) will be read from columne 5 and 6 in the provided TSV file. The resulting fields will be named 'CADD_RawScore' and 'CADD_PhredScore'. If multiple values are provided (comma-separated values), the maximum value will be chosen. Constrained region scores will be read from a BED file, SpliceAI SNP scores will be read from the 'SpliceAI_max' and 'SpliceAI_DP' info fields from the configured VCF file. For SpliceAI_max, the maximum value will be chosen if multiple values are provided. For SpliceAI_DP the values will simply be copied as is. ClinVar IDs, prefixed by the strng 'CV' will be added from the configured VCF file. The cretaed TSV file will contain the configured data columns (and additional standard columns such as position,ref allele, etc.). If multiple values are detected in a comma-separated list then the maximum value will be chosen (except for SpliceAI_SNP_SpliceAI_DP).


Output files
------------

  ::

    test
    |-- annotate_vcf.log # log file
    `-- cohort_mini
        |-- cohort_mini.GRCh38.anno+fil.vcf.gz # annotated VCF file
        |-- cohort_mini.GRCh38.anno+fil.vcf.gz.tbi
        |-- cohort_mini.GRCh38.final.tsv.gz    # TSV file with one entry per variant
    

annotate_cohort
===============

This is a python pipeline for extracting pedigree specific data from annotated VCF files (e.g., created with the annotate_vcf pipeline described above) containing small (SNV, INDEL) and large variants (SVs) and integrating it with various data sources:

* Pedigree data (PED files)
* GADO data 
* Exomiser data
* Additional gene annotations (TSV file)
* GREEN_DB regulatory regions
* Gencode gene annotations

Additional processing steps/processing details:

* Gene symbols will be mapped to current names using a conifgured alias table
* Variants may be filtered for maximum population AF (seperate thresholds for small/large variants)
* Known variants (variants with entries in the ID section of the VCF) are always retained and written to an extra table
* Most severe consequence per variant will be calculated from SnpEff ANN data and configured SO term deleteriousness scores
* Variant locations relative to annotated genes are calculated by overlap with gencode annotations
* Ids of overlapping GREEN_DB regulatory regions are read from Reg_id info fields (e.g., annotated by annotate_vcf)
* Supported inheritance models (recessive, dominant, denovo, any) are calculated from 'inheritance_support' counts:

    - To calculate inheritance_support (=#individuals/GT that support the model-#individuals/GT that contradict the model), the following rules are applied:
        - recessive: 
               -1 for all unaffected sample containing a HOM call,
               +1 for all affected samples with HOM calls that are inherited from mum & dad ,
               0 for all other samples
        - dominant:
               -1 for all unaffected samples with a HOM/HET call,
               -1 for all affected samples containing a HOM call,
               +1 for all affected samples with a HET call that was inherited from an affected sample if available,
               0 for all other samples
        - de_novo:
               -1 for all unaffected samples with a HET/HOM call,
               -1 for all affected samples that inherit a call,
               +1 for all affected samples containing a HET call with high GQ that was not inherited from mum&dad,
               0 for all other samples
        - NOTE that for missing data (including calls with low GQ) we assume that GT that supports the respective inheritance model independent of each other which may lead to the situation that different genotypes for the alleles are assumed per inheritance model.

Example call
------------

  ::

     python cohort_varan/cohort_varan.py \
        annotate_cohort \
        --conf annotate_cohort.example_config.json \
        --threads 1 \
        --out test_cohort/

Config file
-----------
JSON file with the following structure (cf. example below):

* ``dataset_name``: used to create output file names/directories
* ``input_data``: list of considered pedigrees, links to input data files (SNV/CNV VCFs, Pedigree files, GADO+Exomiser files, alias table, gene annotation GFF3 file, GREEN_DB data table, Table of additional gene annotations, gene symbol alias table
* ``filters``: filter thresholds
* ``output_fields``: info fields copied form input VCFs to output tables
* ``d_score_calc``: configuration of how input info fields will be normalized/summarised
* ``so_term`` : SO terms and associated deleteriousness scores and association scores/types. Used for calculation of most severe consequence


Example JSON config file
------------------------

`<annotate_cohort.example_config.json>`_


Output files
------------

  ::

    test_cohort/
    |-- annotate_cohort.log # log file
    |-- cohort_mini.003Neo001.v2r.comphet.tsv.gz    # comphet candidate pairs
    |-- cohort_mini.003Neo001.v2r.genes.tsv.gz      # affected genes and 
    |-- cohort_mini.003Neo001.v2r.known_vars.tsv.gz # known variants
    |-- cohort_mini.003Neo001.v2r.vars.tsv.gz       # annotated/filtered variants
    |-- cohort_mini.v2r.effective_conf.json         # used configuration
    |-- cohort_mini.v2r.idx.tsv.gz                  # index with links to input data for this cohort per pedigree
    |-- cohort_mini.v2r.idx.tsv.gz.tbi
    |-- cohort_mini.v2r.stats.tsv.gz                # Filter statistics per chromosome/variant type
    `-- cohort_mini.v2r.stats.tsv.gz.tbi
    

