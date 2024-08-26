#!/usr/bin/env python3
"""
This script will add pass or fail to the table based on the filters used.
These tags can be used to filter the tabel once more guideliness have been set.
The tags are used to create a log file as well which reports how many primers passed or failed each filter.
The output of this will include the full table a log file and a filtered table.
"""

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description="give arguments to filter script")
parser.add_argument("-i", nargs = 1, required=True, help="file with sequence info example: DMAS-0.txt")
parser.add_argument("-p", nargs = 1, required=True, help="file with all designed primers")
parser.add_argument("-s", nargs = 1, required=True, help="specificity filter: off, strict, loose")
parser.add_argument("-S", nargs=1, required=True, help="SNP filter: off, strict, loose")
parser.add_argument("-t", nargs=1, required=True, help="Secondary structure filter: off, strict, loose")
parser.add_argument("-v", nargs=1, required=True, help="Validation filter off, strict, loose")
parser.add_argument("-l", nargs=1, required=True, help="log file")
args = parser.parse_args()

# get the arguments
input_file = args.i[0]
primers_file = args.p[0]
log_file = args.l[0]
specificity_filter = args.s[0]
SNP_filter = args.S[0]
secondary_structure_filter = args.t[0]
validation_filter = args.v[0]

# Open the table and save the contents
with open(primers_file,'r') as table:
    lines = table.readlines()
table.close()

# get the sequence ID
seq_ID = os.path.splitext(args.i[0])[0]
seq_ID = seq_ID.split("_")[0]
seq_ID = seq_ID.split("/")[-1]
# write a new table
name = seq_ID + "_primers.tsv"
full_table = open(primers_file, 'w')
# create a filtered table
filtered_table_loose = open(seq_ID + "_filtered_loose.tsv", 'w')
filtered_table_strict = open(seq_ID + "_filtered_strict.tsv", 'w')

# go through the table
line_nr = 0
for line in lines:
    line = line.strip().split('\t')
    # Ignore the header and write a new one
    if line_nr == 0:
        full_table.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength_common\tSNPs_FWD\tSec_str_FWD\tSNPs_REV\tSec_str_REV\tDeltaG_template\tFWD_validation\tREV_validation\tAmplicon\tAmp_length\tPredicted_structure\tAmplicon_delta_G\tForward_specificity\tReverse_specificity\tSpecificity_filter_loose\tSpecificity_filter_strict\tSNP_filter_loose\tSNP_filter_strict\tSec_str_filter_loose\tSec_str_filter_strict\tValidation_filter_loose\tValidation_filter_strict\n")
        filtered_table_loose.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength_common\tSNPs_FWD\tSec_str_FWD\tSNPs_REV\tSec_str_REV\tDeltaG_template\tFWD_validation\tREV_validation\tAmplicon\tAmp_length\tPredicted_structure\tAmplicon_delta_G\n")
        filtered_table_strict.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength_common\tSNPs_FWD\tSec_str_FWD\tSNPs_REV\tSec_str_REV\tDeltaG_template\tFWD_validation\tREV_validation\tAmplicon\tAmp_length\tPredicted_structure\tAmplicon_delta_G\n")
        line_nr += 1
        continue
    ################################################################################################
    ###################################   Specificity   ############################################
    ################################################################################################
    # rules from https://academic.oup.com/clinchem/article/59/10/1470/5622018?login=true fig 6 are used as filter criteria
    # The criteria are:
    # Strict: 
    # - no mismatches in either primer = off-target
    # - primers with at least 4 mismatches for a single primer = no off-target
    # - primers with a total of at least 5 mismatches between both primers = no off-target
    # Loose:
    # - primers with at least 3 mismatches for a single primer = no off-target
    # - primers with a total of at least 4 mismatches between both primers = no off-target

    # MM primer1	MM primer2	sum MM	loose			strict
    # 0				0			0		off-target		off-target
    # 1				0			1		off-target		off-target
    # 0				1			1		off-target		off-target
    # 2				0			2		off-target		off-target
    # 0				2			2		off-target		off-target
    # 1				1			2		off-target		off-target
    # 2				1			3		off-target		off-target
    # 1				2			3		off-target		off-target
    # 3				0			3		no off-target	off-target
    # 0				3			3		no off-target	off-target
    # 1				3			4		no off-target	off-target
    # 3				1			4		no off-target	off-target
    # 2				2			4		no off-target	off-target
    # 4				0			4		no off-target	no off-target
    # 0				4			4		no off-target	no off-target
    # 2				3			5		no off-target	no off-target
    # 3				2			5		no off-target	no off-target
    # 4				1			5		no off-target	no off-target
    # 1				4			5		no off-target	no off-target
    # 3				3			6		no off-target	no off-target

    # OFF
    specifity_lost = 0
    Specificity_tag_loose = "off"
    Specificity_tag_strict = "off"

    if specificity_filter != "off":
        try:
            # get the number of mismatches
            Forward_specificity = line[23]
            Forward_specificity = Forward_specificity.replace("[","").replace("]","").replace("'","")
            Forward_specificity = Forward_specificity.split(":")
            Reverse_specificity = line[24]
            Reverse_specificity = Reverse_specificity.replace("[","").replace("]","").replace("'","")
            Reverse_specificity = Reverse_specificity.split(":")
            # get the original input coordinates
            input_file = open(input_file, 'r').readline()
            coord = input_file.split("\t")[-2]
            input_file.close()
        except:
            Specificity_tag_loose = "error: could not get specificity"
            Specificity_tag_strict = "error: could not get specificity"

        # LOOSE
        def specificity_loose(nm_forward, nm_reverse):
            # 3 on a single primer to pass
            if nm_forward < 3 and nm_reverse < 3:
                # total of 4 mismatches to pass
                if nm_forward == 2 and nm_reverse == 2:
                    Specificity_tag_loose = "PASS"
                else:
                    Specificity_tag_loose = "FAIL"
            else:
                Specificity_tag_loose = "PASS"
            return Specificity_tag_loose

        def specificity_strict(nm_forward, nm_reverse):
            # total of at least 5 mismatches to pass
            if nm_forward + nm_reverse < 5:
                # 4 on a single primer to pass
                if (nm_forward == 4 and nm_reverse == 0) or (nm_forward == 0 and nm_reverse == 4):
                    Specificity_tag_strict = "PASS"
                else:
                    Specificity_tag_strict = "FAIL"
            else:
                Specificity_tag_strict = "PASS"
            return Specificity_tag_strict
        
        # placeholder values
        nm_forward_prev = 10
        nm_reverse_prev = 10
        # The find the lowest amount of mismatches in off-targets: if to low fail
        for i in range(0,10):
            try:
                chr_forward = Forward_specificity[i]
                start_forward = Forward_specificity[i+1]
                nm_forward = int(Forward_specificity[i+3])
            except:
                continue

            try:
                chr_reverse = Reverse_specificity[i]
                start_reverse = Reverse_specificity[i+1]
                nm_reverse = int(Reverse_specificity[i+3])
            except:
                continue

            # check if not on target:
            if chr_forward == coord.split(":")[0] and chr_reverse == coord.split(":")[0]:
                if int(start_forward) > int(coord.split(":")[1].split("-")[0]) and int(start_forward) < int(coord.split(":")[1].split("-")[1]) and int(start_reverse) > int(coord.split(":")[1].split("-")[0]) and int(start_reverse) < int(coord.split(":")[1].split("-")[1]):
                    Specificity_tag_loose = "PASS"
                    Specificity_tag_strict = "PASS"
                # On the same chromosome but not the same region
                else:
                    if nm_forward_prev < nm_forward:
                        nm_forward = nm_forward_prev
                    if nm_reverse_prev < nm_reverse:
                        nm_reverse = nm_reverse_prev
                    if Specificity_tag_loose != "FAIL":
                        Specificity_tag_loose = specificity_loose(nm_forward, nm_reverse)
                    if Specificity_tag_strict != "FAIL":
                        Specificity_tag_strict = specificity_strict(nm_forward, nm_reverse)
                        
            # off target
            else:
                if nm_forward_prev < nm_forward:
                    nm_forward = nm_forward_prev
                if nm_reverse_prev < nm_reverse:
                    nm_reverse = nm_reverse_prev
                if Specificity_tag_loose != "FAIL":
                    Specificity_tag_loose = specificity_loose(nm_forward, nm_reverse)
                if Specificity_tag_strict != "FAIL":
                    Specificity_tag_strict = specificity_strict(nm_forward, nm_reverse)


    ################################################################################################
    ###################################   SNPs   ##################################################
    ################################################################################################
    
    # if filter is off
    SNP_tag_loose = "off"
    SNP_tag_strict = "off"
    SNP_lost = 0

    if SNP_filter != "off":
        if SNP_filter != "off":
            # get the found SNPs
            try:
                SNPs_FWD = line[12]
                # none found
                if SNPs_FWD == "0 found":
                    SNPs_FWD = 0
                # get a list of all the positions found
                else:
                    SNPs_FWD_list = SNPs_FWD.replace("{","").replace("}","").split(":")
                    SNPs_FWD = []
                    for i in range(0,len(SNPs_FWD_list),2):
                        SNPs_FWD.append(int(SNPs_FWD_list[i]))
                SNPs_REV = line[14]
                # none found
                if SNPs_REV == "0 found":
                    SNPs_REV = 0
                # get a list of all the positions found
                else:
                    SNPs_REV_list = SNPs_REV.replace("{","").replace("}","").split(":")
                    SNPs_REV = []
                    for i in range(0,len(SNPs_REV_list),2):
                        SNPs_REV.append(int(SNPs_REV_list[i]))
            # failed to get the SNPs
            except:
                SNP_tag_loose = "error: could not get SNPs"
                SNP_tag_strict = "error: could not get SNPs"

            # LOOSE
            def SNP_loose(SNPs, SNP_tag_loose):
                # if the SNPs are not on position -2,-3 or -4 they pass
                for i in SNPs:
                    if i in [-2,-3,-4]:
                        SNP_tag = "FAIL"
                # if it has not yet failed it passes
                if SNP_tag_loose != "FAIL":
                    SNP_tag_loose = "PASS"
                return SNP_tag_loose

            
            # STRICT
            def SNP_strict(SNPs_forward, SNPs_reverse):
                # if there are any SNPs in the primers: fail
                if SNPs_FWD == [] and SNPs_REV == []:
                    SNP_tag_strict = "PASS"
                else:
                    SNP_tag_strict = "FAIL"
                return SNP_tag_strict

            # loop the found SNPs for that line in the table
            if SNPs_FWD == 0 and SNPs_REV == 0:
                SNP_tag_loose = "PASS"
                SNP_tag_strict = "PASS"
            else:
                # strict
                SNP_tag_strict = SNP_strict(SNPs_FWD, SNPs_REV)
                # loose
                if "_F_" in line[0]:
                    SNP_tag_loose = SNP_loose(SNPs_FWD, SNP_tag_loose)
                elif "_R_" in line[0]:
                    SNP_tag_loose = SNP_loose(SNPs_REV, SNP_tag_loose)

        ################################################################################################
        ###################################   Secondary structure   ####################################
        ################################################################################################
        
        # if off
        Sec_str_tag_loose = "off"
        Sec_str_tag_strict = "off"
        Sec_str_lost = 0

        # secondary structure templates
        try:
            # forward positions
            Sec_str_temp_FWD = line[13]
            if Sec_str_temp_FWD == "0 predicted":
                Sec_str_temp_FWD = 0
            else:
                Sec_str_temp_FWD = Sec_str_temp_FWD.replace("[","").replace("]","").split(",")
            # deltaG
            delta_G = float(line[16])
            # reverse positions
            Sec_str_temp_REV = line[15]
            if Sec_str_temp_REV == "0 predicted":
                Sec_str_temp_REV = 0
            else:
                Sec_str_temp_REV = Sec_str_temp_REV.replace("[","").replace("]","").split(",")
            # amplicon
            amplicon = line[21]
            amp_delta_G = float(line[22])
        except:
            Sec_str_tag_loose = "error: could not get secondary structure templates"
            Sec_str_tag_strict


        # check the template

        # if the delta G is less than -15 it fails if any structures are found.
        if delta_G < -15.0:
            if (Sec_str_temp_FWD != 0 and Sec_str_temp_REV != 0) or (Sec_str_temp_FWD !=[] and Sec_str_temp_REV != []):
                Sec_str_tag_strict = "FAIL_template"
            else:
                Sec_str_tag_strict = "PASS"
        else:
            Sec_str_tag_strict = "PASS"
        

        if float(delta_G) < -15:
            for i in Sec_str_temp_FWD:
                if i in [-1,-2,-3,-4,-5,-6]:
                    Sec_str_tag_loose = "FAIL_template"
                else:
                    Sec_str_tag_loose = "PASS"
        else:
            Sec_str_tag_loose = "PASS"


        # check amplicon
        if amp_delta_G < -15.0:
            Sec_str_tag_loose = "FAIL_amplicon"
            Sec_str_tag_strict = "FAIL_amplicon"
        else:
            if amp_delta_G < -5.0:
                if "(" in  amplicon or ")" in amplicon:
                    Sec_str_tag_loose = "FAIL_amplicon"
                    Sec_str_tag_strict = "FAIL_amplicon"
                if Sec_str_tag_loose != "FAIL_template":
                    Sec_str_tag_loose = "PASS"
                if Sec_str_tag_strict != "FAIL_template":
                    Sec_str_tag_strict = "PASS"
            if Sec_str_tag_loose != "FAIL_template":
                Sec_str_tag_loose = "PASS"
            if Sec_str_tag_strict != "FAIL_template":
                Sec_str_tag_strict = "PASS"

    ################################################################################################
    ###################################   Validation   #############################################
    ################################################################################################
    """
    I would not be too stringent on this filter as the validation is not a certainty it depends on a lot of primer3 settings.
    It's more an indication
    """
    # if off
    validation_tag_loose = "off"
    validation_tag_strict = "off"
    validation_lost = 0

    if validation_filter != "off":
        # get the validation results
        try:
            validation_FWD = line[17]
            validation_REV = line[18]
        except:
            validation_tag_loose = "error: could not get validation results"
            validation_tag_strict = "error: could not get validation results"
        
        # strict
        if validation_FWD != "considered 1, ok 1" or validation_REV != "considered 1, ok 1":
            validation_tag_strict = "FAIL"
        else: 
            validation_tag_strict = "PASS"

        
        # loose
        if "_F_" in line[0]:
            if validation_FWD != "considered 1, ok 1":
                validation_tag_loose = "FAIL"
            else:
                validation_tag_loose = "PASS"
        elif "_R_" in line[0]:
            if validation_REV != "considered 1, ok 1":
                validation_tag_loose = "FAIL"
            else:
                validation_tag_loose = "PASS"

    ################################################################################################
    ###################################   Write the table   ########################################
    ################################################################################################

    # Full table
    full_table.write("\t".join(line) + "\t" + Specificity_tag_loose + "\t" + Specificity_tag_strict + "\t" + SNP_tag_loose + "\t" + SNP_tag_strict + "\t" + Sec_str_tag_loose + "\t" + Sec_str_tag_strict + "\t" + validation_tag_loose + "\t" + validation_tag_strict + "\n")
    line_nr += 1

    # drop the last 2 columns
    line = line[:-2]
    #Filtered table
    if Specificity_tag_loose== "PASS" and SNP_tag_loose == "PASS" and Sec_str_tag_loose == "PASS" and validation_tag_loose == "PASS":
        filtered_table_loose.write("\t".join(line) + "\n")
    if Specificity_tag_strict == "PASS" and SNP_tag_strict == "PASS" and Sec_str_tag_strict == "PASS" and validation_tag_strict == "PASS":
        filtered_table_strict.write("\t".join(line) + "\n")

full_table.close()
filtered_table_loose.close()
filtered_table_strict.close()

################################################################################################
# Append to the log file
################################################################################################

# open the log file
log = open(log_file, "a")

# total number of lines = total number of pirmers
total = line_nr - 1

# open the table with pandas
table = pd.read_csv(primers_file, sep="\t")
# Fails on Specificity_filter
specifity_lost_loose = len(table[table["Specificity_filter_loose"] == "FAIL"])
specificity_lost_strict = len(table[table["Specificity_filter_strict"] == "FAIL"])
# Fails on SNP_filter
SNP_lost_loose = len(table[table["SNP_filter_loose"] == "FAIL"])
SNP_lost_strict = len(table[table["SNP_filter_strict"] == "FAIL"])
# Fails on Sec_str_filter
Sec_str_lost_loose = len(table[table["Sec_str_filter_loose"] == "FAIL_template"])
Sec_str_lost_strict = len(table[table["Sec_str_filter_strict"] == "FAIL_template"])
Sec_str_lost_amp_loose = len(table[table["Sec_str_filter_loose"] == "FAIL_amplicon"])
Sec_str_lost_amp_strict = len(table[table["Sec_str_filter_strict"] == "FAIL_amplicon"])
# replace
table["Sec_str_filter_loose"] = table["Sec_str_filter_loose"].replace("FAIL_template", "FAIL")
table["Sec_str_filter_strict"] = table["Sec_str_filter_strict"].replace("FAIL_template", "FAIL")
table["Sec_str_filter_loose"] = table["Sec_str_filter_loose"].replace("FAIL_amplicon", "FAIL")
table["Sec_str_filter_strict"] = table["Sec_str_filter_strict"].replace("FAIL_amplicon", "FAIL")
# Fails on Validation_filter
validation_lost_loose = len(table[table["Validation_filter_loose"] == "FAIL"])
validation_lost_strict = len(table[table["Validation_filter_strict"] == "FAIL"])

# replace off with PASS
table["Specificity_filter_loose"] = table["Specificity_filter_loose"].replace("off", "PASS")
table["Specificity_filter_strict"] = table["Specificity_filter_strict"].replace("off", "PASS")
table["SNP_filter_loose"] = table["SNP_filter_loose"].replace("off", "PASS")
table["SNP_filter_strict"] = table["SNP_filter_strict"].replace("off", "PASS")
table["Sec_str_filter_loose"] = table["Sec_str_filter_loose"].replace("off", "PASS")
table["Sec_str_filter_strict"] = table["Sec_str_filter_strict"].replace("off", "PASS")
table["Validation_filter_loose"] = table["Validation_filter_loose"].replace("off", "PASS")
table["Validation_filter_strict"] = table["Validation_filter_strict"].replace("off", "PASS")

# remaining (all passed)
passed_loose = len(table[(table["Specificity_filter_loose"] == "PASS") & (table["SNP_filter_loose"] == "PASS") & (table["Sec_str_filter_loose"] == "PASS") & (table["Validation_filter_loose"] == "PASS")])
passed_strict = len(table[(table["Specificity_filter_strict"] == "PASS") & (table["SNP_filter_strict"] == "PASS") & (table["Sec_str_filter_strict"] == "PASS") & (table["Validation_filter_strict"] == "PASS")])

# append to the log file
log.write(seq_ID + "\t" + str(total) + "\t" + str(specifity_lost_loose) + "\t" + str(specificity_lost_strict) + "\t" + str(SNP_lost_loose) + "\t" + str(SNP_lost_strict) + "\t" + str(Sec_str_lost_loose) + "\t" + str(Sec_str_lost_strict) + "\t" + str(Sec_str_lost_amp_loose) + "\t" + str(Sec_str_lost_amp_strict) + "\t" + str(validation_lost_loose) + "\t" + str(validation_lost_strict) + "\t"+ str(passed_loose) + "\t" + str(passed_strict) + "\n")
# close log file
log.close()