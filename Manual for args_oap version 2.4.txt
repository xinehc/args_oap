# args_oap_for_conda
Source code of package `args_oap`

## installation
+ via conda (only available for `linux-64, python=3.7`)

```bash
conda install -c bioconda -c xiaole99 args_oap=2.3.7
```
please note that at this moment only python=3.7 is support, if python!=3.7, you may want to create a new conda environment:
    
```bash
conda create -n args_oap -c xiaole99 -c bioconda args_oap=2.3.7 python=3.7
source activate args_oap
```

+ via source

args_oap depends on `diamond>=0.9.24`, `minimap2`, `fastp`, `bbmap`, `samtools`, `blast`. If your system has all the dependencies, then:
```bash
git clone -b test https://github.com/xiaole99/args_oap_for_conda
cd args_oap_for_conda
python setup.py install # use python3 if needed
```

## example
Two toy examples (100k paired-end reads, 100bp each) are provided in `example/inputdir`:

```bash
# git clone -b test https://github.com/xiaole99/args_oap_for_conda
# cd args_oap_for_conda

args_oap stage_one -i example/inputfqs -m example/meta-data.txt -o example/output -f 'fa' -n 8
args_oap stage_two -i example/output/extracted.fa -m example/output/meta_data_online.txt -o example/output -n 8
```

## todo: need to modify the path in the original perl files, otherwise cannot run
Now the perl files use the binaries in the `bin` folder, need to change the path e.g. `./bin/diamond -> diamond` in future updates
