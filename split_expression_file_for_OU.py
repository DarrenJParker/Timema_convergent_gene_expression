# split_expression_file_for_OU.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options needed. To see help specifing -h')


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** split_expression_file_for_OU.py | Written by DJP, 05/10/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("This script takes the logCPM files from 10sp_edgeR.R and produces files that can be used to run Ornstein-Uhlenbeck (OU) models in the R package ouch") 
		print("\n**** USAGE **** \n")

		print("python3.5 split_expression_file_for_OU.py \n\n")
		 

		sys.exit(2)

	else:
		print("i dont know")
		sys.exit(2)


###############################################################################
### infiles

WB_infile_name = "./Gene_exp_for_OU_models/WB_10sp_mean_CPM.csv"
RT_infile_name = "./Gene_exp_for_OU_models/RT_10sp_mean_CPM.csv"
LG_infile_name = "./Gene_exp_for_OU_models/LG_10sp_mean_CPM.csv"

### out dir

WB_GE_for_OU_dir = "WB_GE_for_OU"
RT_GE_for_OU_dir = "RT_GE_for_OU"
LG_GE_for_OU_dir = "LG_GE_for_OU"

try:
	os.mkdir(WB_GE_for_OU_dir)
except OSError as exc:
	print("\n********* WARNING ************\n" + WB_GE_for_OU_dir + " already exists. files may be overwritten.")

try:
	os.mkdir(RT_GE_for_OU_dir)
except OSError as exc:
	print("\n********* WARNING ************\n" + RT_GE_for_OU_dir + " already exists. files may be overwritten.")

try:
	os.mkdir(LG_GE_for_OU_dir)
except OSError as exc:
	print("\n********* WARNING ************\n" + LG_GE_for_OU_dir + " already exists. files may be overwritten.")



###

WB_conv_gene_set = set(["OG-1014", "OG-1039", "OG-1145", "OG-118" , "OG-1246", "OG-1269", "OG-1355", "OG-1487", "OG-1493", "OG-1578",
						"OG-1638", "OG-1651", "OG-1703", "OG-1849", "OG-1891", "OG-1915", "OG-1928", "OG-1945", "OG-1953", "OG-2051",
						"OG-2066", "OG-210" , "OG-2155", "OG-2188", "OG-2197", "OG-2300", "OG-243" , "OG-2440", "OG-2530", "OG-256",
						"OG-260" , "OG-2616", "OG-2632", "OG-2675", "OG-2713", "OG-2829", "OG-2847", "OG-2854", "OG-2856", "OG-2874",
						"OG-2893", "OG-2898", "OG-2928", "OG-2999", "OG-410" , "OG-44"  , "OG-622" , "OG-649" , "OG-663" , "OG-697", 
						"OG-708" , "OG-730" , "OG-763" , "OG-831" , "OG-94"  , "OG-961" , "OG-987"])



RT_conv_gene_set = set(["OG-100" , "OG-1004", "OG-1023", "OG-1026", "OG-1028", "OG-1031", "OG-1073", "OG-1090", "OG-1112", "OG-1116",
						"OG-1126", "OG-1142", "OG-1183", "OG-1191", "OG-1195", "OG-1210", "OG-1211", "OG-1228", "OG-1236", "OG-1238",
						"OG-124" , "OG-1241", "OG-1246", "OG-1291", "OG-1292", "OG-1305", "OG-1323", "OG-134" , "OG-1340", "OG-1353",
						"OG-136" , "OG-1373", "OG-139" , "OG-1394", "OG-1397", "OG-1410", "OG-1419", "OG-1432", "OG-1438", "OG-144",
						"OG-1446", "OG-1448", "OG-1464", "OG-1465", "OG-1478", "OG-148" , "OG-1488", "OG-1506", "OG-1521", "OG-1542",
						"OG-1553", "OG-156" , "OG-1567", "OG-1571", "OG-1610", "OG-162" , "OG-1651", "OG-1655", "OG-1668", "OG-167",
						"OG-1670", "OG-1681", "OG-1684", "OG-1706", "OG-1720", "OG-1749", "OG-1776", "OG-179" , "OG-1794", "OG-1811",
						"OG-1835", "OG-184",  "OG-1841", "OG-1861", "OG-1879", "OG-1888", "OG-1911", "OG-1926", "OG-193" , "OG-1937",
						"OG-1965", "OG-1983", "OG-1991", "OG-1993", "OG-2002", "OG-2021", "OG-2072", "OG-2076", "OG-2082", "OG-2089",
						"OG-2105", "OG-2113", "OG-2155", "OG-2156", "OG-2168", "OG-2176", "OG-2189", "OG-2197", "OG-2216", "OG-2223",
						"OG-2236", "OG-2238", "OG-2243", "OG-2267", "OG-2274", "OG-23"  , "OG-2365", "OG-2368", "OG-2398", "OG-2428",
						"OG-2429", "OG-2441", "OG-2444", "OG-2464", "OG-247" , "OG-2472", "OG-2475", "OG-2488", "OG-2504", "OG-2514",
						"OG-2518", "OG-2524", "OG-2525", "OG-2613", "OG-2658", "OG-2686", "OG-2691", "OG-2719", "OG-2731", "OG-2741",
						"OG-2744", "OG-2753", "OG-2763", "OG-2802", "OG-2808", "OG-283" , "OG-2909", "OG-2931", "OG-2944", "OG-2952",
						"OG-2967", "OG-2992", "OG-314" , "OG-315" , "OG-321" , "OG-328" , "OG-349" , "OG-350" , "OG-354" , "OG-355",
						"OG-36"  , "OG-366" , "OG-374" , "OG-375" , "OG-378" , "OG-386" , "OG-394" , "OG-401" , "OG-406" , "OG-414",
						"OG-421" , "OG-438" , "OG-442" , "OG-444" , "OG-445" , "OG-492" , "OG-507" , "OG-511" , "OG-513" , "OG-518",
						"OG-55"  , "OG-557" , "OG-559" , "OG-565" , "OG-631" , "OG-637" , "OG-643" , "OG-652" , "OG-66"  , "OG-670",
						"OG-689" , "OG-705" , "OG-707" , "OG-712" , "OG-721" , "OG-738" , "OG-757" , "OG-758" , "OG-770" , "OG-775",
						"OG-779" , "OG-805" , "OG-810" , "OG-819" , "OG-89"  , "OG-915" , "OG-930" , "OG-939" , "OG-963" , "OG-964",
						"OG-967" , "OG-976" , "OG-995"])


LG_conv_gene_set = set(["OG-1000", "OG-1007", "OG-1017", "OG-1034", "OG-1073", "OG-1081", "OG-1114", "OG-1128", "OG-1132", "OG-1146",
						"OG-1147", "OG-1148", "OG-1179", "OG-1202", "OG-1208", "OG-1238", "OG-1240", "OG-1247", "OG-1274", "OG-1289",
						"OG-1326", "OG-1337", "OG-134" , "OG-1371", "OG-1377", "OG-1379", "OG-1380", "OG-1385", "OG-1402", "OG-1449",
						"OG-145" , "OG-1464", "OG-1473", "OG-1486", "OG-1524", "OG-1528", "OG-1551", "OG-1592", "OG-1605", "OG-1651",
						"OG-1656", "OG-1660", "OG-1662", "OG-1663", "OG-1664", "OG-1666", "OG-1681", "OG-1684", "OG-1685", "OG-1687",
						"OG-1690", "OG-1695", "OG-1699", "OG-1703", "OG-1707", "OG-1711", "OG-1712", "OG-1729", "OG-1742", "OG-1782",
						"OG-180" , "OG-1818", "OG-1821", "OG-1824", "OG-1839", "OG-1840", "OG-1843", "OG-1852", "OG-1877", "OG-1894",
						"OG-1899", "OG-1942", "OG-1948", "OG-1969", "OG-1990", "OG-2006", "OG-2029", "OG-2031", "OG-2045", "OG-2048",
						"OG-205" , "OG-2080", "OG-2085", "OG-2091", "OG-2099", "OG-2117", "OG-2118", "OG-2122", "OG-2128", "OG-2130",
						"OG-2146", "OG-2155", "OG-2172", "OG-2191", "OG-2197", "OG-2212", "OG-2221", "OG-2232", "OG-2236", "OG-2259",
						"OG-226" , "OG-2260", "OG-2277", "OG-2285", "OG-2288", "OG-2296", "OG-2311", "OG-233" , "OG-2341", "OG-2365",
						"OG-237" , "OG-2376", "OG-2407", "OG-2408", "OG-2417", "OG-2440", "OG-2457", "OG-2474", "OG-2503", "OG-2504",
						"OG-2524", "OG-2537", "OG-2538", "OG-2550", "OG-2599", "OG-2608", "OG-262" , "OG-2623", "OG-2647", "OG-2662",
						"OG-2675", "OG-2676", "OG-2689", "OG-2692", "OG-2702", "OG-2707", "OG-2716", "OG-2728", "OG-2730", "OG-2738",
						"OG-2764", "OG-2768", "OG-2793", "OG-2839", "OG-2840", "OG-2853", "OG-2855", "OG-2861", "OG-2891", "OG-2916",
						"OG-2921", "OG-2924", "OG-2931", "OG-2971", "OG-2973", "OG-2985", "OG-2995", "OG-325" , "OG-341" , "OG-363",
						"OG-376" , "OG-382" , "OG-386" , "OG-406" , "OG-445" , "OG-452" , "OG-470" , "OG-480" , "OG-488" , "OG-501",
						"OG-502" , "OG-518" , "OG-54"  , "OG-559" , "OG-563" , "OG-57"  , "OG-571" , "OG-594" , "OG-617" , "OG-64",
						"OG-642" , "OG-651" , "OG-667" , "OG-714" , "OG-721" , "OG-725" , "OG-745" , "OG-761" , "OG-763" ,  "OG-77",
						"OG-772" , "OG-774" , "OG-796" , "OG-808" , "OG-834" , "OG-84"  , "OG-856" , "OG-872" , "OG-89"  , "OG-903",
						"OG-911" , "OG-937" , "OG-939" , "OG-942",  "OG-944" , "OG-971"])



line_N = 0
WB_infile = open(WB_infile_name)

sample_order = []

for line in WB_infile:

	line = line.rstrip("\n").split(",")
	line_N = line_N + 1

	if line_N == 1:
		for i in line:
			sample_order.append(i.strip('"').split("_")[0])
	else:
		gene_name = line[0].strip('"')
		if gene_name in WB_conv_gene_set:
			outfile = open(os.path.join(WB_GE_for_OU_dir, gene_name  + ".txt"), "w")
			
			for i in range(0, len(line)):
				outfile.write(sample_order[i] + "\t" + line[i].strip('"') + "\n")
				

			outfile.close()	
				
			
line_N = 0
RT_infile = open(RT_infile_name)

sample_order = []

for line in RT_infile:

	line = line.rstrip("\n").split(",")
	line_N = line_N + 1

	if line_N == 1:
		for i in line:
			sample_order.append(i.strip('"').split("_")[0])
	else:
		gene_name = line[0].strip('"')
		if gene_name in RT_conv_gene_set:
			outfile = open(os.path.join(RT_GE_for_OU_dir, gene_name  + ".txt"), "w")
			
			for i in range(0, len(line)):
				outfile.write(sample_order[i] + "\t" + line[i].strip('"') + "\n")
				

			outfile.close()	
				
			
line_N = 0
LG_infile = open(LG_infile_name)

sample_order = []

for line in LG_infile:

	line = line.rstrip("\n").split(",")
	line_N = line_N + 1

	if line_N == 1:
		for i in line:
			sample_order.append(i.strip('"').split("_")[0])
	else:
		gene_name = line[0].strip('"')
		if gene_name in LG_conv_gene_set:
			outfile = open(os.path.join(LG_GE_for_OU_dir, gene_name  + ".txt"), "w")
			
			for i in range(0, len(line)):
				outfile.write(sample_order[i] + "\t" + line[i].strip('"') + "\n")
				

			outfile.close()	
				
			


print("\n\nFinished, Portia\n\n")
	
	














