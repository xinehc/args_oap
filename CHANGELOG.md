## Changelog
#### Version 3.2.4 (27. Oct, 2023)
+ Linting.
+ Add output `blastout.filtered.txt` and `extracted.filtered.fa`.
+ Skip invalid samples in `stage_one` instead of raise exceptions.
+ Fix a bug causing 16S calculation being slow and slightly overestimated.
+ Fix a bug about TPM normalization.

#### Version 3.2.3 (26. June, 2023)
+ Add a gene capO
+ Reduce peak meamory usage for large files https://github.com/xinehc/args_oap/issues/29.

#### Version 3.2.2 (11. January, 2023)
+ Fix naming issue https://github.com/xinehc/args_oap/issues/15.

#### Version 3.2.1 (22. December, 2022)
+ Fix fastq reading bug.
+ Rename `scov` to `copy` in output files.

#### Version 3.2 (16. October, 2022)
+ Simplify interface.
+ Remove dependencies on manually entered `metadata.txt`.
+ Support bioconda installation.

#### Version 3.1.4 (11. October, 2022)
+ Remove a duplicated gene ID (ARCH67_P638154538) in KO30 name list.
+ Remove gene msrE, mphE and some multidrug genes from database.
+ Add RPKM/TPM normalization.

#### Version 3.1.3 (09. September, 2022)
+ Update database to the release version (29082022 short).
+ Add parameter `-s` for skipping 16S/cells calculation in `stageone`.
+ Fix a bug when empty \_2 file is being used as input.

#### Version 3.1.2 (24. August, 2022)
+ Support customized database (testing).
+ Fix a bug of ppm normalization, now the formula is `#hits * 1e6 / #reads`.
+ Simply stagetwo's script, add more information for users.

#### Version 3.1.1 (13. June, 2022)
(If you are expecting results in the unit of ARGs copies per cell (from essential single copy marker genes), this version of fixation won't affect you.)
+ Fix a bug related to the calculation of 16S rDNA copies. Now the numerator is the aligned length instead of read length, and the denominator is the subject length instead of 1432 (the averaged length of 16S rDNA). This will lead to a slight drop of #16S in meta\_data\_online.txt. 
+ Fix a bug when estimating 16S rDNA numbers, which causes some 16S to be counted more than one times when multiple high scoring pairs (HSPs) are returned by blastn.

#### Version 3.1 (09. June, 2022)
+ Minor changes to the SARG database (see [SARG v3.0-M](https://smile.hku.hk/pipeline/#/Indexing/download)).
+ Update diamond to the latest version (from 0.9.24 to 2.0.15), add a new parameter `-w` (query coverage, default 65%) and modify the default parameter of `-v` (identity, default 40%) in stageone to compensate the difference in USCMG estimation.
+ Add a `mt_mode` switcher in stagetwo's blastx to make the pipeline faster when more than 5 cores are available.
+ Add a logger to make stagetwo's stdout more clear.

#### Version 3.0 (04. June, 2022)
+ We updated the SARG database and the corresponding structure file to version 3.0 ([SARG v3.0-M](https://smile.hku.hk/pipeline/#/Indexing/download)) .
+ We dropped bbmap and usearch from the pipeline, now args_oap support both linux and osx.
+ We modified the 16S estimation process by changing minimap2 to bwa + blastn, as minimap2 does not work well for reads that are super short (e.g. below 100 bp, see [https://github.com/lh3/minimap2/issues/363#issuecomment-473387994](https://github.com/lh3/minimap2/issues/363#issuecomment-473387994)).
+ We fixed the version of diamond to 0.9.24 (and python to 3.7.\*), as the latest version of diamond (2.0.15) will give ~10% more hits of USCMGs. The sensitivity of the newer version of diamond is under evaluation. We hope to remove this constraint in future updates.
+ Bug fixed:
    + Fixed a bug that caused the worst hits (instead of the best) to be picked in stagetwo's blastx when multiple candidates of ARGs can be found.
    + Fixed a bug that caused some multi-component ARGs hits to be uncounted in stagetwo's aggregation process.
    + Fixed a bug in stageone that caused USCMG to be slightly overestimated.
    + Fixed a bug in stageone that caused parameters `-x` `-y` `-v` to be ignored.
