# PCN-db-pipeline
## Requirements: biopython, kallisto, SRA-Toolkit, pysradb, and ncbi-datasets-cli

### This program can either be run on your Desktop or Duke Compute Cluster (DCC). DCC is suggested to use due to the large amount of data that is downloaded from the thousands of samples. 

#### Make a top-level directory with three directories inside, named "data", "results", and "src". Now copy all source code files in this repository into "src".

#### The files needed to run the program with are found in the data folder (ncbidataset_assembly and ncbidataset_refseq). The data from these files were obtained by downloading a list of all sequenced baterial genome from https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2&assembly_level=3%3A3. The ref seq and assembly ID columns were isolated. In each file, there are 36,623 samples.

#### The first steps of this pipeline, convert these sample IDs into run accession IDs and create the reference genome path. Two txt files, a table with run accession IDs and a list of the genome paths  will be created. Then, the sample's reference genome will be downloaded if a run accession ID is present, these genomes will be downloaded into the ref-genomes folder in results. There should be 8,612 genomes downloaded, as that is the number of samples that have a run ID.
#### A metadata csv is then created, a master list of all the samples were are downloading data for. Next, fastq files are downloaded for each run ID and put into the SRA folder in results.
#### The next three functions parse through the SRA and reference genomes using kallisto and create files for the three kallisto folders in results.
#### Lastly the final tables of results are created as "chromosome_plasmid_copy_numbers.csv", "ARG_copy_numbers.csv", and "replicon_lengths.csv".

### How to run on DCC:
#### To run on DCC, in your conda environment install bioconda, kallisto, ncbi-datasets-cli and pysradb: 
##### conda install -c bioconda kallisto
##### conda install -c conda-forge ncbi-datasets-cli
##### pip install pysradb
#### You will then need to download the SRA-Toolkit module that is available on DCC:
##### module load SRA-Toolkit

