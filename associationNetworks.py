import os, sys, getopt
import pandas as pd
import numpy as np
import scipy.stats as ss
import time

def usage():
   print 'Usage:'
   print '\tpython associationNetwork.py --mbx path_mbx --mgx path_mgx --Pthreshold 1e-10\n'
   print 'type:'
   print '\tpython associationNetwork.py --help'
   print 'for more details\n'
   sys.exit(2)

def parse_parameter():
   try:
       opts, args = getopt.getopt(sys.argv[1:],shortopts="",longopts=["help","mbx=","mgx=","pThreshold=","output=","","","",""])
   except getopt.GetoptError:
       usage()
   if len(opts) == 0:
       usage()
   hash_config = dict()
   for opt, arg in opts:
       if opt == '--help':
           print 'python associationNetwork.py --mbx path_mbx --mgx path_mgx --pThreshold 1e-10\n'
           print '--mbx\t{ full path to metabolomics table, with features as rows and samples as columns, feature and sample names should be given, tsv file}\n'
           print '--mgx\t{ full path to metagenomics table, with features as rows and samples as columns, feature and sample names should be given, tsv file}\n'
           print '--output\t{ output directory }\n'
           print '--test\t { available options are  fisher, spearman, pearson }\n'
           print '--pThreshold\t{ threshold of significant feature pair }\n'
	   print '--mbx_abundance\t{ abundance threshold for binarization of mbx table to perform Fisher exact test, default 0}\n'
	   print '--mgx_abundance\t{ abundance threshold for binarization of mgx table to perform Fisher exact test, default 0}\n'
           sys.exit(0)
       elif opt == "--mbx":
	   path_mbx = os.path.abspath(arg)
	   if not os.path.isfile(path_mbx):
	       print 'metabolomics file \n' + path_mbx + '\ndoes not exist'
	       sys.exit(2)
	   hash_config['mbx'] = path_mbx
       elif opt == "--mgx":
	   path_mgx = os.path.abspath(arg)
	   if not os.path.isfile(path_mgx):
	       print 'metagenomics file \n' + path_mgx + '\ndoes not exist'
	       sys.exit(2)
           hash_config['mgx'] = path_mgx
       elif opt == "--output":
	   dir_output = os.path.abspath(arg)
	   if os.path.isdir(dir_output):
               print 'output folder \n' + dir_output + '\n already exists'
	       sys.exit(2)
           hash_config['output'] = dir_output
       elif opt == "--test":
	   if arg in ['fisher','spearman','pearson']:
	       hash_config['test'] = arg
           else:
               print 'please give another test name\n'
               sys.exit(2)
       elif opt == "--pThreshold":
           pThreshold = float(arg)
           hash_config['pThreshold'] = pThreshold
       elif opt == "--mbx_abundance":
           mbx_abundance = float(arg)
           hash_config['mbx_abundance'] = mbx_abundance
       elif opt == "--mgx_abundance":
           mgx_abundance = float(arg)
           hash_config['mgx_abundance'] = mgx_abundance
       else:
	   usage()
   return hash_config

def fisher(hash_config):
   fw_out_p = open(hash_config['output'] + 'association_fisher.tsv',"w")
   df_table1 = pd.read_csv(hash_config['mgx'],index_col=0,sep="\t")
   df_table2 = pd.read_csv(hash_config['mbx'],index_col=0,sep="\t")
   list_intersection = list( set(df_table1.columns.values).intersection(set(df_table2.columns.values)) )
   df_table1 = df_table1.loc[:,list_intersection]
   df_table2 = df_table2.loc[:,list_intersection]
   df_table1 = df_table1.ix[df_table1.sum(axis=1) > hash_config['mgx_abundance'],:]
   df_table2 = df_table2.ix[df_table2.sum(axis=1) > hash_config['mbx_abundance'],:]
   fw_out_p.write("\t" + "\t".join(df_table2.index.values) + "\n" )
   if df_table1.shape[0]==0 or df_table2.shape[0]==0:
       print("empty table\n")
       sys.exit(0)
   list_row_r = [0]*df_table2.shape[0]
   list_row_p = [0]*df_table2.shape[0]
   n_col = len(list_intersection)
   print(n_col)

   for index_table1, row_table1 in df_table1.iterrows():
       #tic = time.clock()
       count = 0
       for index_table2, row_table2 in df_table2.iterrows():
           a = ( (row_table1.values + row_table2.values) >= 2 ).sum()
           b = row_table2.values.sum() - a
           c = row_table1.values.sum() - a
           d = n_col - a - b - c
           #print(a,b,c,d)
           list_row_r[count], list_row_p[count] = ss.fisher_exact([[a,b],[c,d]])
           count += 1
       #toc = time.clock()
       #print( toc - tic )
       fw_out_p.write( index_table1 + "\t" + "\t".join(map(str,list_row_p)) + "\n" )
   fw_out_p.close()

def spearman(hash_config):
   fw_out_p = open(hash_config['output'] + 'association_spearman_p.tsv',"w")
   fw_out_r = open(hash_config['output'] + 'association_spearman_r.tsv',"w")
   df_table1 = pd.read_csv(hash_config['mgx'],index_col=0,sep="\t")
   df_table2 = pd.read_csv(hash_config['mbx'],index_col=0,sep="\t")
   list_intersection = list( set(df_table1.columns.values).intersection(set(df_table2.columns.values)) )
   df_table1 = df_table1.loc[:,list_intersection]
   df_table2 = df_table2.loc[:,list_intersection]
   fw_out_p.write("\t" + "\t".join(df_table2.index.values) + "\n" )
   fw_out_r.write("\t" + "\t".join(df_table2.index.values) + "\n" )

   if df_table1.shape[0]==0 or df_table2.shape[0]==0:
       print("empty table\n")
       sys.exit(0)
   list_row_r = [0]*df_table2.shape[0]
   list_row_p = [0]*df_table2.shape[0]
   n_col = len(list_intersection)
   for index_table1, row_table1 in df_table1.iterrows():
       tic = time.clock()
       count = 0
       for index_table2, row_table2 in df_table2.iterrows():
           v1 = row_table1.values + np.random.randn(1,n_col)*1e-13
           v2 = row_table2.values + np.random.randn(1,n_col)*1e-13
           list_row_r[count], list_row_p[count] = ss.mstats.spearmanr(v1,v2)
           count += 1
       toc = time.clock()
       print( toc - tic )
       fw_out_p.write( index_table1 + "\t" + "\t".join(map(str,list_row_p)) + "\n" )
       fw_out_r.write( index_table1 + "\t" + "\t".join(map(str,list_row_r)) + "\n" )
   fw_out_p.close()
   fw_out_r.close()

def pearson(hash_config):
   fw_out_p = open(hash_config['output'] + 'association_pearson_p.tsv',"w")
   fw_out_r = open(hash_config['output'] + 'association_pearson_r.tsv',"w")
   df_table1 = pd.read_csv(hash_config['mgx'],index_col=0,sep="\t")
   df_table2 = pd.read_csv(hash_config['mbx'],index_col=0,sep="\t")
   list_intersection = list( set(df_table1.columns.values).intersection(set(df_table2.columns.values)) )
   df_table1 = df_table1.loc[:,list_intersection]
   df_table2 = df_table2.loc[:,list_intersection]
   fw_out_p.write("\t" + "\t".join(df_table2.index.values) + "\n" )
   fw_out_r.write("\t" + "\t".join(df_table2.index.values) + "\n" )

   if df_table1.shape[0]==0 or df_table2.shape[0]==0:
       print("empty table\n")
       sys.exit(0)
   list_row_r = [0]*df_table2.shape[0]
   list_row_p = [0]*df_table2.shape[0]
   n_col = len(list_intersection)
   for index_table1, row_table1 in df_table1.iterrows():
       tic = time.clock()
       count = 0
       for index_table2, row_table2 in df_table2.iterrows():
           v1 = row_table1.values + np.random.randn(1,n_col)*1e-13
           v2 = row_table2.values + np.random.randn(1,n_col)*1e-13
           list_row_r[count] = np.corrcoef(v1, v2)[0, 1]
           list_row_p[count] = ss.stats.pearson(v1,v2)
           count += 1
       toc = time.clock()
       print( toc - tic )
       fw_out_p.write( index_table1 + "\t" + "\t".join(map(str,list_row_p)) + "\n" )
       fw_out_r.write( index_table1 + "\t" + "\t".join(map(str,list_row_r)) + "\n" )
   fw_out_p.close()
   fw_out_r.close()

def main():
    hash_config = parse_parameter()
    os.system('mkdir ' + dir_output)
    if hash_config['test'] == 'fisher':
       fisher(hash_config)
    elif hash_config['test'] == 'spearman':
       spearman(hash_config)
    elif hash_config['test'] == 'pearson':
       pearson(hash_config)


if __name__ == "__main__":

