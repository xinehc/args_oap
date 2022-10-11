import os
import sys
import argparse
import pandas as pd
import subprocess

from utils import count_bp, logger

############################################ Arguments and declarations ##############################################
workingdir = os.path.abspath(os.path.dirname(__file__))


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="the potential arg reads from stage one",
                    type=str, default='extracted.fa',metavar='extracted.fa')
parser.add_argument("-m",
                    help="meta data online from stage one",
                    type=str, default='meta_data_online.txt',metavar='meta_data_online.txt')
parser.add_argument("-o",
                    help="Output prefix",
                    type=str, default='output_prefix',metavar='output_prefix')
# optional files input and output
parser.add_argument("-n",
                    help="number of threads used for blastx, default 1",
                    type=int, default='4',metavar='4')
parser.add_argument('-l',
                    help="length filtering default 75 percent",
                    metavar="75",
                    default="75", type=float)
parser.add_argument('-a',
                    help="absolute length filtering default 25",
                    metavar="25",
                    default="25", type=float)
parser.add_argument("-e",
                    help="evalue filtering default 1e-7",
                    type=float, default='1e-7',metavar='1e-7')
parser.add_argument('-d',
                    help="identity filtering default 80",
                    metavar="80",
                    default="80", type=float)
parser.add_argument("-b",
                    help="if set then process the blastx results directly [default off], useful if user want to accelerate the stage two by running blastx paralell",
                    type=str, metavar='extracted.out.tab',default="None")
parser.add_argument("-db",
                    help="reference ARG database which is blastxable after makeblastdb, default is SARG database",
                    type=str, default=workingdir+"/DB/SARG.fasta", metavar='SARG.fasta')
parser.add_argument("-struc1",
                    help="reference ARG database structure, default is SARG database",
                    type=str, default=workingdir+"/DB/structure_20210519_for_pipeline.txt", metavar='structure_20210519_for_pipeline.txt')
parser.add_argument("-struc2",
                    help="reference ARG database structure, default is SARG database",
                    type=str, default=workingdir+"/DB/structure_20210519_for_pipeline.txt", metavar='structure_20210519_for_pipeline.txt')
parser.add_argument("-struc3",
                    help="reference ARG database structure, default is SARG database",
                    type=str, default=workingdir+"/DB/structure_20210519_for_pipeline.txt", metavar='structure_20210519_for_pipeline.txt')


################################################### Programme #######################################################

args = parser.parse_args()

print ("\
------------------------------------------------------------------------\n\
This pipeline is designed to process multisamples ARG identification, this is the part two pipeline\n\
perl $0 <Extracted_fasta> <Meta_data_info> <Catergory> <lenth> <e-value> <identity> <institution name> <email address> <Taskname> <PDF16s> <PDFCELL> <TABLE1> <TABLE2> <TABLE3>\n\
Author: Xiaole Yin\n\
Copyright: Prof. Tong Zhang, The University of Hong Kong\n\
Email: yinlele99@gmail.com\n\
Citation: \n\
1. Yin, Xiaole, Xiao-Tao Jiang, Benli Chai, Liguan Li, Ying Yang, James R. Cole,James M. Tiedje, and Tong Zhang.\n\
ARGs-OAP v2.0 with an expanded SARG database and Hidden Markov Models for enhancement characterization \n\
and quantification of antibiotic resistance genes in environmental metagenomes.Bioinformatics 34, no. 13 (2018): 2263-2270.\n\
2. Yang Y, Jiang X, Chai B, Ma L, Li B, Zhang A, Cole JR, Tiedje JM, Zhang T: ARGs-OAP: online analysis\
pipeline for antibiotic resistance genes detection from metagenomic data using an integrated \
structured ARG-database. Bioinformatics 2016. (optional: antibiotic resistance database)\n\
------------------------------------------------------------------------\n\
")

_cols = ["qseqid","sseqid","pident","length","evalue","bitscore","slen","qlen"]
blast6out = str(args.o)+ ".blastout.txt" 
lenmatch = args.l / 100
evaluematch = args.e 
identitymatch =args.d

if os.path.exists(args.db + '.pdb'):
    dbtype = 'prot'
elif os.path.exists(args.db + '.ndb'):
    dbtype = 'nucl'
else:
    logger.critical('Cannot detect the type of DB. Please run <makeblastdb> first!')
    sys.exit(2)

### process meta data
meta = pd.read_csv(args.m,sep="\t")
if 'CellNumber' not in meta or '#of16Sreads' not in meta or '#ofReads' not in meta:
    logger.critical('Wrong metadata file. Please use <meta-data-online.txt> returned by stageone')
    sys.exit(2)
elif any(meta['CellNumber'] == 0) or any(meta['#of16Sreads'] == 0) or any(meta['#ofReads'] == 0):
    logger.critical('Some sample has empty #reads/#16s/#cells, cannot normalize!')
    sys.exit(2)

### blast
if args.b != 'None':
    blast6out = args.b
else:
    nbp, nseq = count_bp(args.i, ">")
    mt_mode = "1" if nbp / int(args.n) >= 2500000 else "0"

    logger.info("Extracting target sequences ({} reads, {} threads, mode {}) ...".format(nseq, args.n, mt_mode))
    blast = 'blastx' if dbtype == 'prot' else 'blastn'
    cmd = [
        blast,
        "-db", args.db, 
        "-query", args.i, 
        "-out", blast6out, 
        "-outfmt", " ".join(["6"]+_cols),
        "-max_target_seqs", "5",
        "-evalue", str(evaluematch),
        "-num_threads", str(args.n), 
        "-mt_mode", mt_mode
        ]
    subprocess.run(cmd, check=True)

### process blast6out
logger.info("Merging files ...")
df = pd.read_csv(blast6out, sep="\t", header=None, names=_cols)

if len(df) == 0:
    logger.critical('No target sequence detected!')
    sys.exit(2)

df = df[( df["pident"] >= identitymatch)  & ( df["evalue"] <= evaluematch)]

## count max qlen
perc = (df['qlen'] == df['qlen'].max()).sum()/len(df) * 100
lenmatchabs = args.a if dbtype =='prot' else args.a * 3
if perc != 100:
    logger.warning('Reads have uneven lengths, please double check.')

logger.info('Read length: {} ({}%), absolute alignment length cutoff in aa: {}'.format(
    df['qlen'].max(), perc, args.a))

df = df[df['length'] >= lenmatchabs]
df["scov"] = df["length"] / df["slen"] 
df["qcov"] = df["length"] * 3 / df["qlen"] if dbtype == 'prot' else df["length"] / df["qlen"]
df = df[df["qcov"] >= lenmatch]

df = df.sort_values(["qseqid", "evalue", "bitscore", "length"], ascending=[True, True, False, False])
df = df.drop_duplicates(subset='qseqid', keep="first")
df['Name'] = df['qseqid'].str.rsplit('_', 1).str.get(0)
df['sseqid'] = df['sseqid'].astype(str)

### process SARG structure
if 'DB/single-component_structure.txt' in args.struc1:
    struc_list = []
    for file, count in zip([args.struc1, args.struc2, args.struc3], [1, 0.333, 0.5]):
        struc_list.append(pd.read_csv(file, sep="\t", dtype = str).assign(count = count))
    struc = pd.concat(struc_list, ignore_index=True)
else:
    struc = pd.read_csv(args.struc1, sep='\t', dtype = str).assign(count = 1)

levels = struc.columns[:-1]

### merge blast6out results with meta and SARG structure#
if not all(x in set(struc[levels[0]]) for x in df["sseqid"].unique()):
    logger.warning("Not all extracted sequences can be found in the structure file, please check the database and the structure file.")


result = pd.merge(df, struc, left_on = 'sseqid', right_on = levels[0], how = 'inner')
result['scov'] = result['scov'] * result['count'] # account for different weight

## calculate length normalized count
result['rpk'] = result['count']/ (result['slen'] / 1000)

logger.info("Saving output files ...")
for level in levels:
    for measure, normalizer, name in zip(['scov', 'scov', 'count', 'rpk', 'rpk'], 
                                         ['CellNumber', '#of16Sreads', '#ofReads', '#ofReads', '#ofReads'],
                                         ['normalize_cellnumber', 'normalize_16s', 'ppm', 'fpkm', 'percentage']):
        agg = result.groupby([level, 'Name'])[measure].sum().reset_index()
        out = pd.merge(agg, meta, on = 'Name', how ='outer')
        out['value'] = out[measure]/out[normalizer]

        if name in {'ppm', 'fpkm'}:
            out['value'] = out['value'] * 1e6

        if name in {'percentage'}:
            out = pd.merge(out, out.groupby('Name')['value'].sum().reset_index(), on='Name')
            out['value'] = out['value_x'] / out['value_y']

        out.set_index(level).pivot(columns = 'Name', values = 'value').fillna(0).sort_index().to_csv(
            str(args.o) + '.' + name + '.' + level + '.txt', sep ='\t') 

        ## save unnormalized subject coverage or count
        if name in {'normalize_16s', 'ppm'}:
            out.set_index(level).pivot(columns = 'Name', values = measure).fillna(0).sort_index().to_csv(
                str(args.o) + '.unnormalize_' + measure + '.' + level + '.txt', sep ='\t') 

logger.info('Done')
