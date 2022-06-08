import os
import argparse
import pandas as pd

from datetime import datetime

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
# parser.add_argument("-db",
#                     help="reference ARG database which is blastxable after makeblastdb, default is SARG database",
#                     type=str, default=workingdir+"/DB/SARG.3.fasta", metavar='SARG.3.fasta')
parser.add_argument("-db",
                    help="reference ARG database which is blastxable after makeblastdb, default is SARG database",
                    type=str, default=workingdir+"/DB/SARG_20210519.fasta", metavar='SARG_20210519.fasta')
parser.add_argument("-fa",
                    help="reference ARG database in fasta format, default is SARG database",
                    type=str, default=workingdir+"/DB/SARG_20210519.fasta", metavar='SARG_20210519.fasta')
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

blast6out = str(args.o)+ ".blastout.txt" 
lenmatch = args.l / 100
evaluematch = args.e 
identitymatch =args.d
begin = datetime.now().strftime("%H:%M:%S")
if args.b != 'None':
    blast6out = args.b
else:
    os.system("blastx -query "+str(args.i)+" -out "+blast6out+" -db "+str(args.db)+" -evalue "+str(evaluematch)+" -num_threads "+str(args.n)+' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -max_target_seqs 5')

### process blast6out #######
df = pd.read_csv(blast6out,sep="\t",header=None,names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","slen","qlen"])
df = df[( df["pident"] >= identitymatch)  & ( df["evalue"] <= evaluematch)]
df["ratio"] = df["length"] / df["slen"] 
df["qratio"] = df["length"] * 3 / df["qlen"] 
df = df[df["qratio"] >= lenmatch]
df = df.sort_values(["qseqid","evalue","bitscore"], ascending=[True, True, False])
df = df.drop_duplicates(subset='qseqid', keep="first")
df['qseqid_modi'] = df['qseqid'].str.rsplit('_', 1).str.get(0)


### process meta data #######
meta = pd.read_csv(args.m,sep="\t")

### process SARG structure ########
struc1 = pd.read_csv(args.struc1,sep="\t")
struc2 = pd.read_csv(args.struc2,sep="\t")
struc3 = pd.read_csv(args.struc3,sep="\t")


### merge blast6out results with meta and SARG structure#########
result1 = pd.merge(df, struc1, left_on='sseqid',right_on='SARG.Seq.ID',how='inner')
result1_type = result1.groupby(['qseqid_modi', 'Type'])[['ratio']].agg('sum')
result1_subtype = result1.groupby(['qseqid_modi', 'Subtype'])[['ratio']].agg('sum')
result1_gene = result1.groupby(['qseqid_modi', 'sseqid'])[['ratio']].agg('sum')
result1["count"]=1
ppm1_type = result1.groupby(['qseqid_modi', 'Type'])[['count']].agg('sum')
ppm1_subtype = result1.groupby(['qseqid_modi', 'Subtype'])[['count']].agg('sum')
ppm1_gene = result1.groupby(['qseqid_modi', 'sseqid'])[['count']].agg('sum')

result2 = pd.merge(df, struc2, left_on='sseqid',right_on='SARG.Seq.ID',how='inner')
result2['ratio_modifi']=result2['ratio']/3
result2_type = result2.groupby(['qseqid_modi', 'Type'])[['ratio_modifi']].agg('sum')
result2_subtype = result2.groupby(['qseqid_modi', 'Subtype'])[['ratio_modifi']].agg('sum')
result2_gene = result2.groupby(['qseqid_modi', 'sseqid'])[['ratio']].agg('sum')
result2["count"]=0.333
ppm2_type = result2.groupby(['qseqid_modi', 'Type'])[['count']].agg('sum')
ppm2_subtype = result2.groupby(['qseqid_modi', 'Subtype'])[['count']].agg('sum')
ppm2_gene = result2.groupby(['qseqid_modi', 'sseqid'])[['count']].agg('sum')

result3 = pd.merge(df, struc3, left_on='sseqid',right_on='SARG.Seq.ID',how='inner')
result3['ratio_modifi']=result3['ratio']/2
result3_type = result3.groupby(['qseqid_modi', 'Type'])[['ratio_modifi']].agg('sum')
result3_subtype = result3.groupby(['qseqid_modi', 'Subtype'])[['ratio_modifi']].agg('sum')
result3_gene = result3.groupby(['qseqid_modi', 'sseqid'])[['ratio']].agg('sum')
result3["count"]=0.5
ppm3_type = result3.groupby(['qseqid_modi', 'Type'])[['count']].agg('sum')
ppm3_subtype = result3.groupby(['qseqid_modi', 'Subtype'])[['count']].agg('sum')
ppm3_gene = result3.groupby(['qseqid_modi', 'sseqid'])[['count']].agg('sum')

result2_type = result2_type.rename(columns={'ratio_modifi':'ratio'})
result2_subtype = result2_subtype.rename(columns={'ratio_modifi':'ratio'})
result2_gene = result2_gene.rename(columns={'ratio_modifi':'ratio'})
result3_type = result3_type.rename(columns={'ratio_modifi':'ratio'})
result3_subtype = result3_subtype.rename(columns={'ratio_modifi':'ratio'})
result3_gene = result3_gene.rename(columns={'ratio_modifi':'ratio'})


result1_type = result1_type.add(result2_type,fill_value=0)
result1_type = result1_type.add(result3_type,fill_value=0).reset_index()
result1_subtype = result1_subtype.add(result2_subtype,fill_value=0)
result1_subtype = result1_subtype.add(result3_subtype,fill_value=0).reset_index()
result1_gene = result1_gene.add(result2_gene,fill_value=0)
result1_gene = result1_gene.add(result3_gene,fill_value=0).reset_index()

ppm1_type = ppm1_type.add(ppm2_type,fill_value=0)
ppm1_type = ppm1_type.add(ppm3_type,fill_value=0).reset_index()
ppm1_subtype = ppm1_subtype.add(ppm2_subtype,fill_value=0)
ppm1_subtype = ppm1_subtype.add(ppm3_subtype,fill_value=0).reset_index()
ppm1_gene = ppm1_gene.add(ppm2_gene,fill_value=0)
ppm1_gene = ppm1_gene.add(ppm3_gene,fill_value=0).reset_index()

result1_type = pd.merge(result1_type,meta,left_on="qseqid_modi",right_on="Name",how="outer")
result1_subtype = pd.merge(result1_subtype,meta,left_on="qseqid_modi",right_on="Name",how="outer")
result1_gene = pd.merge(result1_gene,meta,left_on="qseqid_modi",right_on="Name",how="outer")

result1_type["percell"] = result1_type["ratio"] / result1_type["CellNumber"]
result1_type["per16s"] = result1_type["ratio"] / result1_type["#of16Sreads"]
result1_subtype["percell"] = result1_subtype["ratio"] / result1_subtype["CellNumber"]
result1_subtype["per16s"] = result1_subtype["ratio"] / result1_subtype["#of16Sreads"]
result1_gene["percell"] = result1_gene["ratio"] / result1_gene["CellNumber"]
result1_gene["per16s"] = result1_gene["ratio"] / result1_gene["#of16Sreads"]

ppm1_type = pd.merge(ppm1_type,meta,left_on="qseqid_modi",right_on="Name",how="outer")
ppm1_subtype = pd.merge(ppm1_subtype,meta,left_on="qseqid_modi",right_on="Name",how="outer")
ppm1_gene = pd.merge(ppm1_gene,meta,left_on="qseqid_modi",right_on="Name",how="outer")
ppm1_type["ppm"] = ppm1_type["count"] / ppm1_type["#ofReads"]
ppm1_subtype["ppm"] = ppm1_subtype["count"] / ppm1_subtype["#ofReads"]
ppm1_gene["ppm"] = ppm1_gene["count"] / ppm1_gene["#ofReads"]

typecell = result1_type.pivot(index='Type', columns='qseqid_modi')['percell'].fillna(0).sort_index()
type16s = result1_type.pivot(index='Type', columns='qseqid_modi')['per16s'].fillna(0).sort_index()
typeppm = ppm1_type.pivot(index='Type', columns='qseqid_modi')['ppm'].fillna(0).sort_index()
subtypecell = result1_subtype.pivot(index='Subtype', columns='qseqid_modi')['percell'].fillna(0).sort_index()
subtype16s = result1_subtype.pivot(index='Subtype', columns='qseqid_modi')['per16s'].fillna(0).sort_index()
subtypeppm = ppm1_subtype.pivot(index='Subtype', columns='qseqid_modi')['ppm'].fillna(0).sort_index()
genecell = result1_gene.pivot(index='sseqid', columns='qseqid_modi')['percell'].fillna(0).sort_index()
gene16s = result1_gene.pivot(index='sseqid', columns='qseqid_modi')['per16s'].fillna(0).sort_index()
geneppm = ppm1_gene.pivot(index='sseqid', columns='qseqid_modi')['ppm'].fillna(0).sort_index()


#### write out ###
typecell.to_csv(str(args.o)+".normalize_cellnumber.type.txt",sep='\t')
type16s.to_csv(str(args.o)+".normalize_16s.type.txt",sep='\t')
typeppm.to_csv(str(args.o)+".ppm.type.txt",sep='\t')
subtypecell.to_csv(str(args.o)+".normalize_cellnumber.subtype.txt",sep='\t')
subtype16s.to_csv(str(args.o)+".normalize_16s.subtype.txt",sep='\t')
subtypeppm.to_csv(str(args.o)+".ppm.subtype.txt",sep='\t')
genecell.to_csv(str(args.o)+".normalize_cellnumber.gene.txt",sep='\t')
gene16s.to_csv(str(args.o)+".normalize_16s.gene.txt",sep='\t')
geneppm.to_csv(str(args.o)+".ppm.gene.txt",sep='\t')
