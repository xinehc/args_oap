# ARGs_OAP (beta)
This repository was created by Xiaole Yin (_xiaole99_) and is currently maintained by Xi Chen (_xinhec_). The goal is to make args_oap faster, and easier to run. 

If you have any questions, please contact Xiaole Yin ([yinlele99@gmail.com](yinlele99@gmail.com)).

## Installation
Conda (osx-64/linux-64):
```bash
conda install -c bioconda -c conda-forge xinehc::args_oap=3.1.1
```

We'd suggest to create a new conda environment (here use `-n args_oap` as an example) to avoid potential conflicts of dependencies:
```bash
conda create -n args_oap -c bioconda -c conda-forge xinehc::args_oap=3.1.1
conda activate args_oap
```

Args_oap depends on `python>=3.7`, `diamond>=2.0.15`, `bwa>=0.7.17`, `blast>=2.12`, `samtools>=1.15`, `fastp>=0.23.2`, `pandas`. If your OS has all the dependencies, then it can be built from source:
```bash
git clone https://github.com/xinehc/args_oap.git
cd args_oap
python setup.py install # use python3 if needed
```

## Example
Two examples (100k paired-end reads, 100 bp each) can be found [here](https://dl.dropboxusercontent.com/s/054ufvfahchfk7f/example.tar.gz). The zipped file can be downloaded manually or using `wget`:

```bash
# conda install wget
wget https://dl.dropboxusercontent.com/s/054ufvfahchfk7f/example.tar.gz
tar -xvf example.tar.gz
cd example

# conda activate args_oap
args_oap stage_one -i inputfqs -m meta-data.txt -o output -f fa -n 8
args_oap stage_two -i output/extracted.fa -m output/meta_data_online.txt -o output/output -n 8
```

After `stage_one`, a `meta_data_online.txt` file can be found in `output`. It summarizes the 16s and cell numbers of each samples, for example:

| SampleID | Name     | Category | ReadLength | #ofReads | #of16Sreads      | CellNumber       |
|----------|----------|----------|------------|----------|------------------|------------------|
| 1        | STAS     | ST       | 100        | 200000   | 8.22929787979405 | 2.82910169634418 |
| 2        | SWHAS104 | SWH      | 100        | 200000   | 7.00954780712517 | 3.16376328499407 |

After `stagetwo`, the normalized ARGs copies per 16s/cells or hits/reads will be shown in several `*_normalized_*.txt` files. For example, `output.normalize_16s.type` means:
+ **normalized_16s** - normalized against 16s rRNA copies
+ **type** - Type of ARGs (the hierarchy in SARG is type -> subtype -> gene)

| Type                                | STAS                  | SWHAS104             |
|-------------------------------------|-----------------------|----------------------|
| aminoglycoside                      | 0.016202273302162954  | 0.07023487211759523  |
| bacitracin                          | 0.014243756749154245  | 0.029022199416685355 |
| beta_lactam                         | 0.0                   | 0.07435959439262277  |
| macrolide-lincosamide-streptogramin | 0.0                   | 0.011963015532727397 |
| multidrug                           | 0.0049211431544144806 | 0.015912491818927718 |
| mupirocin                           | 0.002966724847808158  | 0.004557467877130285 |
| quinolone                           | 0.14468645528642882   | 0.043998731935285834 |
| sulfonamide                         | 0.01345207192984009   | 0.06808763199694169  |
| tetracycline                        | 0.004659396079993181  | 0.04969500656817938  |

###  (mandatory) Prepare the meta-data.txt file
(We hope to remove the manual meta-data.txt preparation step in future updates)

To run the stage one pipeline, you need to 
1. Put all your paired-end fastq/fasta files into one directory (notice that the name of your fastq files should be Name_1.fq and Name_2.fq).  
2. Prepare relative meta-data.txt file.  

Tips:     
* You need keep the first and second column's name as SampleID and Name  
* The SampleID are required to be unique numbers counting from 1 to 2 to 3 etc.  
* Category is the classification of your samples into groups and we will colored your samples in PcoA by this information (online stagetwo version only, see below)
* The meta-data table should be separated by tabular for each of the items   
* The Name of each sample should be the fastq file names for your pair-end Illumina sequencing data, your fastq files will automatically be recognized by Name_1.fq and Name_2.fq, so you need to keep the name consistent with your fq file name. (if you files are end with .fastq or .fasta, you need to change them to end with .fq or .fa)  
   
**Please make sure the meta-data file is pure txt format, if you edit the file under windows, using notepad++ and check the end of each line by cliking View-> Show Symbol -> Show All Characters. If the line is end up with CRLF, please remove the CR by replace \r to nothing in the replace dialogue frame.**

The meta-data.txt shall look like this:

SampleID | Name | Category |ReadLength  
---------|------|-------|----  
 1       | STAS | ST       |100  
 2       | SWHAS104 | SWH  |100  

### (optional) For very big data to run stage two  
For users have very big data and prefer complex running:  
1. users run locally by themselves to get the blastx outfmt 6 format resutls by alighment against SARG2.2.  
**A typical scene is that users can paralelly run the blastx on clusters by multi-nodes, and then merge the blastx output as the input for the -b option.**  
2. use -b option for the stage two script:   

```bash
args_oap stage_two -i extracted.fa -m meta_data_online.txt -o output -b merge_blastx.out.txt  
```

### (optional) Stage two pipeline on Galaxy system and download results  
Go to http://smile.hku.hk/SARGs  and using the module ARG_OAP.    
  
1. Using **ARG_OAP** -> **Upload Files** module to upload the extracted fasta file and meta_data_online.txt file generated in stage one into Galaxy    
2. Click **ARG_OAP** and **Ublast_stagetwo**, select your uploaded files    
3. For \"Column in Metadata:\" chose the column you want to classify your samples (default: 3)  
  
Click **Execute** and you can find four output files for your information  
  
After a while or so, you will notice that their are four files generated for your information.    
   
**File 1 and 2**: PcoA figures of your samples and other environment samples generated by ARGs abundance matrix normalization to 16s reads number and cell number    
**File 3 and 4**: Other tabular mother tables which including the profile of ARGs type and sub type information, as long as with other environment samples mother table. File3 results of ARGs abundance normalization against 16S reads number; File 4 results of ARGs abundance normalization against cell number  
  
There are some questions raised by users, please refer to the [FAQ](https://github.com/biofuture/Ublastx_stageone/wiki/FAQ) for details. To run ARG OAP locally, users should download the source code into local computer system (Unix/Linux). Users can upload the generated files for stage two onto our Galaxy analysis platform (http://smile.hku.hk/SARGs) or use the local version of stage two script.   


## Changes log
#### Version 3.1.1 (13. June, 2022)
+ If you are expecting results in the unit of copy of ARGs per cell (from essential single copy marker genes), this version of fixation won't affect you. In detail:
   + A bug related to the calculation of estimated 16S rDNA numbers was fixed. Now the numerator is the aligned length instead of read length, and the denominator is the subject length instead of 1432 (the averaged length of 16S rDNA). This will lead to a slight drop of #16s in meta_data_online.txt. 
   + Fix a bug when estimating cell numbers from the 16S copyrighter. Bugs caused some 16s to be counted more than one times when multiple high scoring pairs (HSPs) are returned by blastn.


#### Version 3.1 (09. June, 2022)
+ Minor changes to the SARG database (see [SARG v3.0-M](https://smile.hku.hk/pipeline/#/Indexing/download)).
+ Update diamond to the latest version (from 0.9.24 to 2.0.15), add a new parameter `-w` (query coverage, default 65%) and modify the default parameter of `-v` (identity, default 40%) in stageone to compensate the difference in USCMG estimation.
+ Add a `mt_mode` switcher in stagetwo's blastx to make the pipeline faster when more than 5 cores are available.
+ Add a logger to make stagetwo's stdout more clear.

#### Version 3.0 (04. June, 2022)
+ We updated the SARG database and the corresponding structure file to version 3.0 ([SARG v3.0-M](https://smile.hku.hk/pipeline/#/Indexing/download)) .
+ We dropped bbmap and usearch from the pipeline, now args_oap support both linux and osx.
+ We modified the 16s estimation process by changing minimap2 to bwa + blastn, as minimap2 does not work well for reads that are super short (e.g. below 100 bp, see [https://github.com/lh3/minimap2/issues/363#issuecomment-473387994](https://github.com/lh3/minimap2/issues/363#issuecomment-473387994)).
+ We fixed the version of diamond to 0.9.24 (and python to 3.7.\*), as the latest version of diamond (2.0.15) will give ~10% more hits of USCMGs. The sensitivity of the newer version of diamond is under evaluation. We hope to remove this constraint in future updates.
+ Bug fixed:
    + Fixed a bug that caused the worst hits (instead of the best) to be picked in stagetwo's blastx when multiple candidates of ARGs can be found.
    + Fixed a bug that caused some multi-component ARGs hits to be uncounted in stagetwo's aggregation process.
    + Fixed a bug in stageone that caused USCMG to be slightly overestimated.
    + Fixed a bug in stageone that caused parameters `-x` `-y` `-v` to be ignored.

---    
**Notice:**  
  
This tools only provide the required scripts for ARGs-OAP1.0/2.0 pipeline  
  
This pipeline is distributed in the hope to achieve the aim of management of antibiotic resistant genes in envrionment, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.This pipeline is only allowed to be used for non-commercial and academic purpose.  
  
**The SARG database is distributed only freely used for academic prupose, any commercial use should require the agreement from the developer team.**   
