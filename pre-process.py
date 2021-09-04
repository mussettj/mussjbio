#!/usr/bin/env

import os
import sys
import collections
import pandas
import re
from datetime import date
from Bio import SeqIO

# passed to the pre-process script.
# `./pre-process.py $file $min_maf $abs_min $min_call_rate $meta $startweek $endweek`;

# RUN USING:
# python pre-process.py shortened_cog_alignment.fasta 0 0 0 shortened_cog_metadata.csv 55 61
# with example shortened cog metadata files.


# sys/argv[0] is the name of the script e.g. pre-process.py

# input COG UK or other alignment file
file = sys.argv[1]
if not os.path.isfile(file):
    raise FileNotFoundError("Alignment file not found")

f = open('%s.pipeline.log' % file, 'a')  # open file and append to it. Creates file if it doesn't exist yet.

# Set this variable to be the minimum minor allele frequency at a position considered real.
# E.g. 0.01 means min of 1/100 reads
min_maf = int(sys.argv[2])
if min_maf == 0:      # if zero, set to minimum maf of 0.01
    min_maf = 0.01

# Need to have an absolute minimum as well, otherwise singletons will be viewed as potential SNPs in small datasets.
abs_min = int(sys.argv[3])
if abs_min == 0:   # if zero, set to minimum of 4
    abs_min = 4

# Set this variable to the minimum proportion of bases at a position with a CATG call (as opposed to ambiguous
# or missing).
min_call_rate = int(sys.argv[4])
if min_call_rate == 0:  # if zero, set to minimum call rate of 0.9
    min_call_rate = 0.9

meta = sys.argv[5]   # cog metadata file needed
if not os.path.isfile(meta):
    raise FileNotFoundError("Metadata file not found")

startweek = int(sys.argv[6])  # start week to search for in pandemic weeks
endweek = int(sys.argv[7])    # end week to search for in pandemic weeks

# PRIORITY VARIANTS TEXT FILE CHECK
# Check for priority variants and make sure these are used for potential panel members if so

if os.path.isfile("priority_variants_pos.txt"):  # if priority variants file exists, read it
    priority = pandas.read_table('priority_variants_pos.txt', header=0)   # read in priority variants file
    prior = []
    banned = []
    priority_variants = []

    priority['subs'] = priority['subs'].apply(str)   # create table from the subs and variant columns
    priority['Variant'] = priority['Variant'].apply(str)
    sub = priority['subs']

    for i in range(len(priority)):  # loop through the table, if the subs column is in the format of Base (C,A,
        # T or G) number and Base, then it is a match
        matches = re.match('([CATG])(\d+)([CATG])', priority.iloc[i, 2])
        if matches:                   # if the subs column matches the regex format, the position is the number and
            # the snp is the second base
            pos = matches.group(2)
            snp = matches.group(3)

            if re.match('(\D)', priority.iloc[i, 0]) and re.match('(\d)', pos) and re.match('([A-Z])', snp) and \
                    priority.iloc[i, 3] > 0:    # if the weight of the snp in the priority file is more than zero,
                # add it to the priority variants list
                prior.append(pos)
            elif priority.iloc[i, 3] < 0:  # if the weight of the snp in the priority file is less than zero, add it
                # to the list of banned snps
                banned.append(pos)

        priority_variants.append(priority.iloc[i, 0])  # concatenate the position and snp to the priority variants list


# METADATA SECTION
metadata = pandas.read_csv(meta, header=0, usecols=['sequence_name', 'epi_week'])  # read in sequence name and epi
# week from metadata file
metadata = metadata.set_index('sequence_name')  # set index to sequence name

sam_count = metadata.groupby(["epi_week"]).size().reset_index(name='Sample count')  # count number of samples per week
sam_count.columns = ['Pandemic week', 'Sample count']
sam_count = sam_count.set_index('Pandemic week')

sam_sum = sam_count['Sample count'].sum()   # count total number of samples in metadata

with open('%s.pipeline.log' % file, 'a') as f:  # write out to log the sample count and number of samples loaded
    f.write(f" \n date ran: {date.today()} \n Metadata file: {meta}  \n {sam_count}. \n Loaded week data for {sam_sum} samples."
            f"\n pre-process_covid_sequences.pl seqfile: {sys.argv[1]} minMAF: {min_maf} ABS-min: {abs_min} "
            f"min call rate: {min_call_rate} metadata: {meta} startweek: {startweek} endweek: {endweek} \n")

# Biopython
record_dict = SeqIO.index(file, "fasta")  # read in the alignment fasta file

listed = list(record_dict)  # list of ids/records from the dictionary of sequences from the fasta

ids = []
good = 0
new_list = []
new_x = []
for x in listed:  # for sample in the alignment file
    if x in metadata.index:   # if the sample is in the metadata file
        week = int(metadata.loc[x, 'epi_week'])     # find the week which it was recorded
        sequence = (record_dict[x].seq)  # the sequence of the associated id/name
        split_seq = list(sequence)   # split sequence of bases into columns in to a list
        if (startweek <= week) & (week <= endweek):  # if week sample taken between or equal to week range
            new_x.append(x)   # add the id to list
            ids.append(x)    # add id to list (to be kept as list)
            good += 1  # add one to number of tallied samples used
            new_list.append(list(sequence))  # add listed sequence to list


print(f"\n {good} sequences from {file} fell between weeks {startweek} and {endweek}\n ")

# Add the priority variants files
unique = set(priority_variants)  # remove duplicate entries of priority variant names from list of priority variants

add = []
pv_count = 0
for i in unique:
    if os.path.isfile("%s.fa" % i):   # if priority variant name is a fasta file in the directory
        rec = "%s.fa" % i  # find record of file
        add = SeqIO.read(rec, "fasta") # read in fasta file using SeqIO
        ids.append(add.id)  # add the sample/variant id
        new_x.append(add.id)
        rec_seq = add.seq   # add the fasta sequence of the sample/variant
        new_list.append(list(rec_seq))  # add fasta sequence to list
        pv_count += 1   # priority variant count add 1

# Now create dataframe of the selected sequences plus the sequences from the priority variants files.

sequence_table = pandas.DataFrame(index=new_x, data=new_list)  # create dataframe of all of the accumulated sample
# sequences and their ids, including priority variants/samples
sequence_table.columns += 1  # counting starts from 0, reset column names of bases to start from base 1 at col 1

types = len(list(ids))   # types is number of samples/variants loaded
minii = types * min_maf  # multiply number of loaded samples by the minimum minor allele frequency to find minimum
# number of bases which can be used to discriminate a SNP from reading errors and random mutations.
if minii < abs_min:
    minii = abs_min

with open('%s.pipeline.log' % file, 'a') as f:
    f.write(f"\n {good} sequences from {file} fell between weeks {startweek} and {endweek}\n ")  # print to log
    f.write(f"\n {pv_count} priority variants added to analysis. \n")
    f.write(f"\n Loaded {types} sequences")

with open('%s_covid_variants.tdt' % file, 'w') as g:   # create covid variants tdt file to be used in processing script
    for each_id in ids:
        g.write(f"\t {each_id}")  # write in the header with each sample id used

    g.write(f"\n")

sequence_table.replace(to_replace='[^CATG]', value='N', regex=True)    # replace any non-base values with N.
# sequence_table.replace(to_replace="-", value='N')
# ensure every position in the dataframe is a CATG, or if not, replace with an N.
# For every line/file in the sequence table,

size = sequence_table.shape  # gauge size of the dataframe
final_output = []
snp_count = 0
pri_snp_count = 0
for i in range(1, size[1]):  # for the length of the first row entry, i.e. all columnns
    bad = 0
    base_count = []
    col = sequence_table[i]  # the base at the column position
    for j in range(0,col.size):     # for each row within the column
        if (col[j] == "N"):  # if the base at that row in the column is N, add 1 to bad
            bad = bad+1  # bad is just the count of Ns?

        base_count.append(col[j])   # for each row within the column, add the base to the lise base count
        base_count_join = ''.join(base_count)  # remove the spaces

    counted = collections.Counter(base_count)   # tally the number of counts per base for each column
    new_counted = counted.most_common(3)  # gives the 2 most common occurances and their counts in order.
    #therefore the value in position zero is the Major allele and the value in position one is the minor allele

    if len(new_counted) == 1:  # if only one base counted, not a snp, continue to next column
        continue

    if len(new_counted) == 2:   # if length of counts is 2
        major_base = new_counted[0][0]   # major base is the most common count
        major_count = new_counted[0][1]

        minor_base = new_counted[1][0]  # minor base is second most common count
        minor_count = new_counted[1][1]

        if major_base == "-":  # ignore - and N's in count, wont be snps
            continue

        if major_base == "N":
            continue

        if minor_base == "-":   # if minor - or N, ignore
            continue

        if minor_base == "N":
            continue

    if len(new_counted) == 3:  # if minor base is N then use the next most common base for count
        major_base = new_counted[0][0]
        major_count = new_counted[0][1]

        minor_base = new_counted[1][0]
        minor_count = new_counted[1][1]

        if major_base == "-":
            continue

        if major_base == "N":
            continue

        if minor_base == "-":
            continue

        if minor_base == "N":
            minor_base = new_counted[2][0]
            minor_count = new_counted[2][1]
            if minor_base == "-":
                continue

    # if($n > $min && 1-($bad/$types) >$min_call_rate && $banned{$i} <1)
    if (minor_count > minii) and (1-(bad/types) > min_call_rate) and (str(i) not in banned):  # if minor count is
        # above the threshold of the minimum minor allele frequency, and is not banned, it is a snp
        snp_count += 1
        with open('%s.covid_variants.tdt' % file, 'a') as g:  # write snp out to output file for process script
            g.write("Base " + str(i) + " " + str(major_base) + str(major_count) + str(minor_base) + str(minor_count) + " " + str(base_count_join) + "\n")

        print("Base " + str(i) + " " + str(major_base) + str(major_count) + str(minor_base) + str( minor_count) + " " + str(base_count_join))

    elif str(i) in prior:   # if the snp position is a priority, then add it to the output for processing
        snp_count += 1
        pri_snp_count += 1

        with open('%s_covid_variants.tdt' % file, 'a') as g:  # write out priority snps to output file
            g.write("Base " + str(i) + " " + str(major_base) + str(major_count) + str(minor_base) + str(
                minor_count) + " " + str(base_count_join) + "\n")

print(f"\n Found {snp_count} SNPs, including {pri_snp_count} priority SNPs. The minor allele frequency used to discriminate the SNPs was {minii}.")
# print how many snps were found