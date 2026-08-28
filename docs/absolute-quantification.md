## Example of Absolute Quantification

Two example FASTA datasets, each containing 100,000 paired-end reads of 100 bp, are available [here](https://dl.dropboxusercontent.com/s/pqgftlo24rfc2rd/example.tar.gz).

The compressed file can be downloaded manually or using `wget`:

```bash
# conda install wget
wget https://dl.dropboxusercontent.com/s/pqgftlo24rfc2rd/example.tar.gz
tar -xvf example.tar.gz
cd example
```

The `example/input` directory contains four files:

```text
STAS_1.fa
STAS_2.fa
SWHAS104_1.fa
SWHAS104_2.fa
```

Here, `STAS` and `SWHAS104` are the sample names. The suffixes `_1` and `_2` indicate the forward and reverse reads, respectively, for paired-end sequencing data.

### Prepare the metadata file

Create a tab-delimited metadata file named `meta.txt` with the following format:

```text
sample	SampleVolume_mL	ElutionVolume_uL	DNAConcentration_ng/uL
STAS	50	60	50
SWHAS104	50	60	50
```

The first line must contain the following header exactly:

```text
sample	SampleVolume_mL	ElutionVolume_uL	DNAConcentration_ng/uL
```

The columns are defined as follows:

- `sample`: Sample name. It must match the sample name used in the input FASTA filenames. In this example, the samples are `STAS` and `SWHAS104`.
- `SampleVolume_mL`: Volume of the original sample used for DNA extraction, in mL.
- `ElutionVolume_uL`: Final elution volume of the extracted DNA, in µL. For example, the referenced protocol uses 60 µL, but this value should be provided by the user because it may vary among extraction kits and protocols.
- `DNAConcentration_ng/uL`: Concentration of the extracted DNA, in ng/µL.

### Run absolute quantification

For absolute quantification, place the metadata file at `input/meta.txt`, or specify its location using the `-m FILE` option.

Run:
```bash
# conda activate args_oap
args_oap stage_one -i input -o output -f fa -t 8 -m input/meta.txt
args_oap stage_two -i output -t 8
```
stage_two does not require the -m argument. Absolute quantification will be calculated automatically if the -m argument was provided in stage_one.

Samples with missing input files, missing metadata, or invalid/non-positive metadata values will be skipped.

After `stage_two`, two additional types of output are generated:

```text
metadata_with_absolute_quantification.txt
absolute_quantification.{type,subtype,gene}.txt
```

The absolute-quantification files report ARG abundance as **total ARG copies per mL of the original sample**.

All existing ARGs-OAP outputs remain unchanged.
