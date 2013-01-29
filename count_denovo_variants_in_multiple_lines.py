#!/usr/bin/python2.6
#
# Finds denovo variants in multiple lines via 2 rounds:
# in the first round potential candidates are collected from .vcf files (samtools pileups)
# in the second round the candidates are evaluated by accessing the .bam files directly
#
# In order to separate the variants present in the mother (-> in several malines) from the unique ones (-> in a single maline)
# only positions that have a variant in one line and solid reference bases in several other lines are accepted 
#
# The difference between SNPs and Indels is that snps are int he 2nd round:
# snps are checked against a raw pileup and a bam without a similar snp base is collected as ref [0,1] or with a similar base as alt [1,0]
# indel are checked against a called vcf and a bam with an indel call is collected as alt [1,0] 
# in the end snps are checked for [1, >10], indels only for [1, ?]

import os
import re
import sys
sys.path.append("/home/andreasw/pymodules")
import multiprocessing as mp
import datetime
import optparse
import samtools_parser
import commands


#####################################################################################################################################################
##########################################################            Options               #######################################################
#####################################################################################################################################################

p = optparse.OptionParser()

# mandatory:
p.add_option("-t", help="the title identifier of the target lines", default= None, dest="line_titles")
p.add_option("--target", help="specify target type of variant: hom/het/indel", default=None, dest="target_type")
p.add_option("--bams", help="location of the bams", default=None, dest="bamsfolder")


# optional:
p.add_option("--parent", help="bam of parental line to be included in the 2nd run", default=None, dest="parent_bam")
p.add_option("--chunkfile", help="vcf of parental line to be included in the 2nd run", default=None, dest="chunk_file")
p.add_option("--debug", help="only run on contig0 for debugging", action = "store_true", default=False, dest="debug")
p.add_option("--celegans", help="run on c.elegans data", action = "store_true", default=False, dest="celegans")
p.add_option("--printmode", help="print intermediate data for debugging", action = "store_true", default=False, dest="printmode")
p.add_option("--rejected", help="print all rejected candidate lines", action = "store_true", default=False, dest="rejected")

opts, args = p.parse_args()

line_titles = opts.line_titles
parent_bam = opts.parent_bam
target_type = opts.target_type
chunk_file = opts.chunk_file
bamsfolder = opts.bamsfolder
debug = opts.debug
celegans = opts.celegans
printmode = opts.printmode
rejected = opts.rejected

print opts

if line_titles == None:
    sys.stderr.write("\nPlease supply at least the target line titles with '-t' (e.g. -t SA), the target type [het/hom/indel] with --target and the bams folder with --bams.\n \
					 e.g. -t MA --target indel -- bams ./bams \
					 \n Options:\n --parent [parent line]\n --chunkfile [chunkfile in (contig start stop) format]\n --debug \n")
    raise SystemExit(1)


##################################################
# threshold values for finding candidate snps

cand_min_dp = 5
cand_max_dp = 30
cand_min_qual = 10 
cand_max_fq = -40
cand_min_fq = 40 # for het search

########################################################################
# threshold values for validating candidate snp positions in other lines

eval_min_dp = 3

#####################################################################################################################################################
##########################################################            Functions               #######################################################
#####################################################################################################################################################

#############################################################################################################
# retrieve a raw base view from a given .bam

def retrieve_raw_base_view_from_bam(conpos, bam, bamsfolder):
	
	con = conpos.split("\t")[0]
	pos = conpos.split("\t")[1]
	
	sam_cmd = "samtools mpileup -r %s:%s-%s -AB -Duf ~/bioinfo/qtl_data_preparation_qc/Pristionchus_Hybrid_assembly.fa %s%s | bcftools view -" % (con, pos, pos, bamsfolder,bam)
	if celegans: sam_cmd = "samtools mpileup -r %s:%s-%s -AB -Duf /home/andreasw/bioinfo/malines/c_elegans/c_elegans.WS170.genomic.fa %s%s | bcftools view -" % (con, pos, pos, bamsfolder,bam) 
	if printmode: print sam_cmd
	
	result = commands.getoutput(sam_cmd)
	fields = result.split("\n")
    
	for row in fields:
	    if row[0] == "C":
		return row

#############################################################################################################
# retrieve a raw base view from a given .bam

def retrieve_all_raw_bases_from_bam(conpos, bam, bamsfolder):
	
	con = conpos.split("\t")[0]
	pos = conpos.split("\t")[1]
	
	sam_cmd = "samtools mpileup -r %s:%s-%s -AB -Df ~/bioinfo/qtl_data_preparation_qc/Pristionchus_Hybrid_assembly.fa %s%s " % (con, pos, pos, bamsfolder,bam)
	if celegans: sam_cmd = "samtools mpileup -r %s:%s-%s -AB -Df /home/andreasw/bioinfo/malines/c_elegans/c_elegans.WS170.genomic.fa %s%s " % (con, pos, pos, bamsfolder,bam) 
	#if printmode: print sam_cmd
	
	result = commands.getoutput(sam_cmd)
	fields = result.split("\n")
    
	for row in fields:
	    if row[0] == "C":
		return row

#############################################################################################################
# retrieve a vcf pileup with snp calling from a given .bam

def retrieve_vcf_from_bam(conpos, bam, bamsfolder):
	
	con = conpos.split("\t")[0]
	pos = conpos.split("\t")[1]
	
	sam_cmd = "samtools mpileup -r %s:%s-%s -AB -Duf ~/bioinfo/qtl_data_preparation_qc/Pristionchus_Hybrid_assembly.fa %s%s | bcftools view -v -" % (con, pos, pos, bamsfolder,bam)
	if celegans: sam_cmd = "samtools mpileup -r %s:%s-%s -AB -Duf /home/andreasw/bioinfo/malines/c_elegans/c_elegans.WS170.genomic.fa %s%s | bcftools view -" % (con, pos, pos, bamsfolder,bam)
	if printmode: print sam_cmd
	
	results = commands.getoutput(sam_cmd)
	fields = results.split("\n")
    
	
	result_rows = [] # the problem is that indels lead to 2 rows being returned with the same conpos, therefore BOTH have to be collected and evaluated
	
	for row in fields:
	    if row[0] == "C":
			result_rows.append(row)
	
	#print result_rows
	
	if len(result_rows) == 1:
		return result_rows[0]
	
	elif len(result_rows) > 1:
		for row in result_rows:
			if "INDEL" in row:
				return row
	else:
		return None

#############################################################################################################
# 1st round: map function to do a focal run first through the vcf variant file to find all good candidates for all lines

def map_find_candidate_snps(focal_filename):
		
	cand_snp_rows = {} # collects the candidate snps in the target MAline
	cand_snp_malines = {} # remembers the MAline that a conpos belonged to

	maline_name = focal_filename.split(".")[0]
	
	with open(focal_filename) as focal_maline_file:		
		for row in focal_maline_file:
			
			if row[0] != "#": 
				var = samtools_parser.VCFrow(row)
				
				if var.is_valid_row():
					
					if debug is False or (debug and not celegans and var.con == 5) or (debug and celegans and var.contig_name == "CHROMOSOME_I"): # only one Contig in debugging mode									
						if not celegans:
							if var.is_ref_known() and var.min_qual(cand_min_qual) and var.min_coverage(cand_min_dp) and var.max_coverage(cand_max_dp): # row is accepted as possible candidate
								if (target_type == "hom" and var.is_hom()) or (target_type == "het" and var.is_het()) or (target_type == "indel" and var.is_indel() and var.max_fq(cand_max_fq)): # correct type of variant
									
									conpos = var.get_conpos()
									cand_snp_rows[conpos] = row
									cand_snp_malines[conpos] = maline_name
									
								
						else:
							if var.is_ref_known() and var.min_qual(cand_min_qual) and var.min_coverage(cand_min_dp) and var.max_coverage(cand_max_dp): # row is accepted as possible candidate
							
								conpos = var.contig_name + "\t" + str(var.get_pos())
								candidate_bam = maline_name.split(".")[0] + ".bam"
								celegans_thresholds_ok = False
									
								if target_type != "indel":
									rawbaserow = retrieve_all_raw_bases_from_bam(conpos, candidate_bam, bamsfolder)
									
									if rawbaserow is not None:
										raw_var = samtools_parser.RawBaseRow(rawbaserow)
										if raw_var.is_valid_celegans_denver_snp(var.get_alt()):
											celegans_thresholds_ok = True
											
									if printmode: print "candidate:", celegans_thresholds_ok, row,
										
								if celegans_thresholds_ok or (target_type == "indel" and  var.is_indel()):
									cand_snp_rows[conpos] = row
									cand_snp_malines[conpos] = maline_name
											
	return [cand_snp_rows, cand_snp_malines, focal_filename]
	
#############################################################################################################
# 1st round: reduce function to reduce candidate snps to the unique ones

def reduce_candidate_snps(list_of_dicts):
	
	cand_snp_rows = list_of_dicts[0]
	cand_snp_malines = list_of_dicts[1]
	
	print "1st run done for", list_of_dicts[2], "found", len(cand_snp_rows), "candidates" 
	
	for conpos in cand_snp_rows:
		if printmode: print cand_snp_rows[conpos],
		if conpos not in all_cand_snp_rows:
			all_cand_snp_rows[conpos] = cand_snp_rows[conpos]
			all_cand_malines[conpos] = 	cand_snp_malines[conpos]
			
		else:
			conpos_to_delete.append(conpos)


#############################################################################################################
# 2nd round: map function for snps to collect data on the candidate snps by accessing the positions via "samtools -r"

def map_evaluate_candidate_snps_in_bams(focal_filename):
		
	cand_snp_comparison = {} # collects the conpos with a ref base in this line

	maline_name = focal_filename.split(".")[0]
		
	for candidate_snp_row in all_cand_snp_rows:
		
		cand_var = samtools_parser.VCFrow(all_cand_snp_rows[candidate_snp_row])
		
		target_row = retrieve_raw_base_view_from_bam(candidate_snp_row, focal_filename, bamsfolder)
		
		if printmode: print "target_row:", target_row
		
		if target_row != None:
			target_var = samtools_parser.RawBCFrow(target_row)

			if target_var.is_ref():
				cand_snp_comparison[cand_var.get_conpos()] = [0,1]
				
			elif cand_var.alt in target_var.get_altbase():
				cand_snp_comparison[cand_var.get_conpos()] = [1,0]

	print "done", focal_filename
	if printmode: print "cand_snp_comparison",cand_snp_comparison
	return cand_snp_comparison

#############################################################################################################
# 2nd round: map function for indels to collect data on the candidate snps by accessing the positions via "samtools -r"

def map_evaluate_candidate_indels_in_bams(focal_filename):
		
	cand_snp_comparison = {} # collects the conpos with a ref base in this line

	maline_name = focal_filename.split(".")[0]
		
	for candidate_snp_row in all_cand_snp_rows:
		
		cand_var = samtools_parser.VCFrow(all_cand_snp_rows[candidate_snp_row])
		target_row = retrieve_vcf_from_bam(candidate_snp_row, focal_filename, bamsfolder)
		#print target_row
		
		if target_row != None:
			target_var = samtools_parser.VCFrow(target_row)

			if target_var.is_indel():
				cand_snp_comparison[cand_var.get_conpos()] = [1,0]

	print "done", focal_filename
	if printmode: print "cand_snp_comparison",cand_snp_comparison
	return cand_snp_comparison

#############################################################################################################
# 2nd round: function to test wether a position lies within a valid chunk

def test_chunk_membership(conpos):
	
	in_chunk = False
	
	conpos = conpos.lstrip("Contig")
	con = int(conpos.split("\t")[0])
	pos = int(conpos.split("\t")[1])
	
	if con in chunks:
	
		contig_chunks = chunks[con]
				
		for interval in contig_chunks:
			if interval[0] <= pos and interval[1] >= pos:
				in_chunk = True
				break
		
	return in_chunk
		
#####################################################################################################################################################
##########################################################   1st round: collect candidates   ########################################################
#####################################################################################################################################################

#############################################################################################################
# setup variables

conpos_to_delete = [] # collects the conpos of candidate snps that have occured more than once in the 1st round
all_cand_snp_rows = {}
all_cand_malines = {}
all_cand_snp_comparisons = {}

cand_snp_comparison_container = []

if chunk_file != None:
	
	chunks = {}
	
	with open(chunk_file) as chunk_input:
		for row in chunk_input:
			row = row.rstrip("\n")
			fields = row.split("\t")
			
			contig = int(fields[0])
			start = int(fields[1])
			stop = int(fields[2])
			
			interval = [start, stop]
			
			if contig in chunks:
				chunks[contig].append(interval)
				
			else:
				chunks[contig] = [interval]
				
			
#############################################################################################################
# collect the vcf files to be used

variant_files = []

for focal_filename in os.listdir(os.curdir):
	if focal_filename.startswith(line_titles) and focal_filename.endswith(".vcf") and "allpos" not in focal_filename and "lines" not in focal_filename and "rawbase" not in focal_filename:
		variant_files.append(focal_filename)

print variant_files

#############################################################################################################
# multiprocessing of 1st round

pool = mp.Pool()

for focal_filename in variant_files:
	#reduce_candidate_snps(map_find_candidate_snps(focal_filename))
	pool.apply_async(map_find_candidate_snps, args = (focal_filename,), callback = reduce_candidate_snps)

pool.close()
pool.join()

#############################################################################################################
# delete non-unique snps

print len(all_cand_snp_rows.keys()), "candidates left before deletion of duplicates"

for conpos in conpos_to_delete:
	if conpos in all_cand_snp_rows:
		if rejected: print all_cand_snp_rows[conpos],
		del all_cand_snp_rows[conpos]
		del all_cand_malines[conpos]

print len(all_cand_snp_rows.keys()), "unique candidates left after duplicate deletion"

for key in all_cand_snp_rows.keys():
	if printmode: print all_cand_malines[key], all_cand_snp_rows[key],
	pass

#############################################################################################################
# delete snps in chunks

if chunk_file is not None:
	for conpos in all_cand_snp_rows.keys():
		if test_chunk_membership(conpos) is False:
			if rejected: print all_cand_snp_rows[conpos]
			del all_cand_snp_rows[conpos]
			del all_cand_malines[conpos]
	
	print len(all_cand_snp_rows.keys()), "unique candidates remain in the valid chunks"

if rejected: 1/0

#####################################################################################################################################################
##########################################################   2nd round: validate candidates   #######################################################
#####################################################################################################################################################

#############################################################################################################
# collect the bamfiles to be used

bams_files = []

if bamsfolder == None:
	bamsfolder = os.curdir

for focal_filename in os.listdir(bamsfolder):
	if focal_filename.startswith(line_titles) and focal_filename.endswith(".bam"):
		bams_files.append(focal_filename)
	
if parent_bam != None:
	bams_files.append(parent_bam)

print bams_files

#############################################################################################################
# multiprocessing of 2nd round

pool = mp.Pool()

if target_type != "indel":
	result = pool.map_async(map_evaluate_candidate_snps_in_bams, bams_files)
else:
	result = pool.map_async(map_evaluate_candidate_indels_in_bams, bams_files)
"""

for bam in bams_files:
	print "starting", bam
	print map_evaluate_candidate_indels_in_bams(bam)
"""
pool.close()
pool.join()
print "2nd round done"

cand_snp_comparison_container = result.get()
#print "cand_snp_comparison_container:", cand_snp_comparison_container

#############################################################################################################
# process stored results

for cand_snp_comparison in cand_snp_comparison_container:
	for conpos in cand_snp_comparison:
		if conpos not in all_cand_snp_comparisons:
			all_cand_snp_comparisons[conpos] = cand_snp_comparison[conpos]
		else:
			all_cand_snp_comparisons[conpos][0] += cand_snp_comparison[conpos][0]
			all_cand_snp_comparisons[conpos][1] += cand_snp_comparison[conpos][1]

#print "all_cand_snp_comparisons:", all_cand_snp_comparisons
	
#############################################################################################################
# print valid snps

if target_type == "het":
	outname = line_titles + "lines_denovo_het_snps" + datetime.datetime.now().strftime("%Y%m%d_%H%M") + ".vcf"
elif target_type == "indel":
	outname = line_titles + "lines_denovo_indels" + datetime.datetime.now().strftime("%Y%m%d_%H%M") + ".vcf"
elif target_type == "hom":
	outname = line_titles + "lines_denovo_hom_snps" + datetime.datetime.now().strftime("%Y%m%d_%H%M") + ".vcf"

output = open(outname, "w")
valid_snp_count = 0

for conpos in all_cand_snp_comparisons:
	
	accepted_variant = False
	
	if all_cand_snp_comparisons[conpos][0] == 1: # its a 1 because each candidate will be found again when the same line will be tested in the 2nd round
		if target_type == "indel" or (target_type != "indel" and all_cand_snp_comparisons[conpos][1] > 10) or celegans:
			accepted_variant = True
			
	
	if 	accepted_variant == True:
		
		target_line = all_cand_malines[conpos]
		target_row = all_cand_snp_rows[conpos]

		valid_snp_count += 1
		
		#output.write(str(all_cand_snp_comparisons[conpos][0]) +"\t"+ str(all_cand_snp_comparisons[conpos][1]) +"\t"+ target_line +"\t"+ target_row + "\n")
		output.write(target_line +"\t"+ target_row + "\n")

		print (target_line +"\t"+ target_row), 

print valid_snp_count, "candidates left after round 2"
