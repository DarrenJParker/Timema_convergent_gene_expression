### OU_summarise_and_FDR.py

import sys
import os
import getopt
import decimal
import statsmodels.stats.multitest as smm

try:
	opts, args = getopt.getopt(sys.argv[1:], 'w:r:l:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


WB_dir_name = None
RT_dir_name = None
LG_dir_name = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** OU_summarise_and_FDR.py | Written by DJP, 05/10/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("\n**** USAGE **** \n")
		print("python3.5 OU_summarise_and_FDR.py -w [whole-body OU output dir] -r [Reproductive tract OU output dir] -l [Leg OU output dir]\n\n") 
			
		sys.exit(2)
		
	elif opt in ('-w'):
		WB_dir_name  = arg
	elif opt in ('-r'):
		RT_dir_name = arg
	elif opt in ('-l'):
		LG_dir_name = arg
	else:
		print("i dont know")
		sys.exit(2)
		
		
### get files and join

WB_pvals = []
WB_gene_ord = []
WB_gene_info = {}

header_txt = ""

path = WB_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(".csv"):
			#print (os.path.join(path, name))
			curr_file = open(os.path.join(path, name))
			line_N = 0
			for line in curr_file:
				line_N = line_N + 1
				if line_N == 1:
					header_txt = line.rstrip("\n")
				if line_N > 1:
					line = line.rstrip("\n").split(",")
					genename = line[0].replace(".", "-")
					g_info = line[1:]
					
					WB_pvals.append(float(line[-1]))
					WB_gene_ord.append(genename)
					WB_gene_info[genename] = g_info

RT_pvals = []
RT_gene_ord = []
RT_gene_info = {}

path = RT_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(".csv"):
			#print (os.path.join(path, name))
			curr_file = open(os.path.join(path, name))
			line_N = 0
			for line in curr_file:
				line_N = line_N + 1
				if line_N > 1:
					line = line.rstrip("\n").split(",")
					genename = line[0].replace(".", "-")
					g_info = line[1:]
					
					RT_pvals.append(float(line[-1]))
					RT_gene_ord.append(genename)
					RT_gene_info[genename] = g_info



LG_pvals = []
LG_gene_ord = []
LG_gene_info = {}

path = LG_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(".csv"):
			#print (os.path.join(path, name))
			curr_file = open(os.path.join(path, name))
			line_N = 0
			for line in curr_file:
				line_N = line_N + 1
				if line_N > 1:
					line = line.rstrip("\n").split(",")
					genename = line[0].replace(".", "-")
					g_info = line[1:]
					
					LG_pvals.append(float(line[-1]))
					LG_gene_ord.append(genename)
					LG_gene_info[genename] = g_info


##### adjust pvals

rej, pval_corr = smm.multipletests(WB_pvals, method="fdr_bh")[:2]
WB_pval_corr_list = list(pval_corr)

rej, pval_corr = smm.multipletests(RT_pvals, method="fdr_bh")[:2]
RT_pval_corr_list = list(pval_corr)

rej, pval_corr = smm.multipletests(LG_pvals, method="fdr_bh")[:2]
LG_pval_corr_list = list(pval_corr)

### output p-adjusted files

header_out = header_txt + ",LRT_OU_brown_FDR"

WB_out = open(WB_dir_name.replace("/", "") + "_wFDR.csv", "w")
WB_out.write(header_out + "\n")

for i in range(0,len(WB_gene_ord)):
	gene_name = WB_gene_ord[i]
	gene_FDR = WB_pval_corr_list[i]
	gene_info = WB_gene_info.get(gene_name)
	gene_info_out = ""
	for el in gene_info:
		gene_info_out = gene_info_out + "," + el
	
	WB_out.write(gene_name + gene_info_out + "," + str(gene_FDR) + "\n")


RT_out = open(RT_dir_name.replace("/", "") + "_wFDR.csv", "w")
RT_out.write(header_out + "\n")

for i in range(0,len(RT_gene_ord)):
	gene_name = RT_gene_ord[i]
	gene_FDR = RT_pval_corr_list[i]
	gene_info = RT_gene_info.get(gene_name)
	gene_info_out = ""
	for el in gene_info:
		gene_info_out = gene_info_out + "," + el
	
	RT_out.write(gene_name + gene_info_out + "," + str(gene_FDR) + "\n")


LG_out = open(LG_dir_name.replace("/", "") + "_wFDR.csv", "w")
LG_out.write(header_out + "\n")

for i in range(0,len(LG_gene_ord)):
	gene_name = LG_gene_ord[i]
	gene_FDR = LG_pval_corr_list[i]
	gene_info = LG_gene_info.get(gene_name)
	gene_info_out = ""
	for el in gene_info:
		gene_info_out = gene_info_out + "," + el
	
	LG_out.write(gene_name + gene_info_out + "," + str(gene_FDR) + "\n")
	


##### make summary table

WB_sig = 0
WB_total = 0
for el in WB_pval_corr_list:
	WB_total = WB_total + 1
	if el < 0.05:
		WB_sig = WB_sig + 1

RT_sig = 0
RT_total = 0
for el in RT_pval_corr_list:
	RT_total = RT_total + 1
	if el < 0.05:
		RT_sig = RT_sig + 1
		
LG_sig = 0
LG_total = 0
for el in LG_pval_corr_list:
	LG_total = LG_total + 1
	if el < 0.05:
		LG_sig = LG_sig + 1
		

print("\nNumber of genes where the OU model fits better than Brownian model (FDR < 0.05):")
print("WB\t" + str(WB_sig) + " (" + str((WB_sig / WB_total * 100)) + " %)")
print("RT\t" + str(RT_sig) + " (" + str((RT_sig / RT_total * 100)) + " %)")
print("LG\t" + str(LG_sig) + " (" + str((LG_sig / LG_total * 100)) + " %)\n\n")

print("Done, Avrana Kern\n\n")




