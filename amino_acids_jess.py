#!/usr/bin/env

import os
import sys
import collections
import pandas
import re
from Bio import SeqIO

# ./flatten_fasta.pl is used to flatten mafft alignment before pre-processing

# passed to the pre-process script.
# `./amino_acids_jess.py $file $min_maf $abs_min $min_call_rate $meta $startweek $endweek`;

# RUN USING:
# python amino_acids_jess.py mafft_500_wuhanref_output 0 0 0 x 0 0
# example with mafft amino acid alignment files.

# sys/argv[0] is the name of the script e.g. pre-process.py

# input mafft alignment file.
file = sys.argv[1]
if not os.path.isfile(file):
    raise FileNotFoundError("Alignment file not found")

f = open('%s.pipeline.log' % file, 'a')  # open file and append to it.
# creates file if it doesn't exist yet.

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
if min_call_rate == 0:    # if zero, set to minimum call rate of 0.9
    min_call_rate = 0.9

meta = sys.argv[5]
if not os.path.isfile(meta):    # metadata file unnecessary, but easier to copy csv script into folder to run script.
    meta = 'cog_2020-09-03_metadata.csv'

startweek = int(sys.argv[6])  # set to zero for aa
endweek = int(sys.argv[7])    # set to zero for aa

# Biopython
# read in file
filedata = None
with open(file, 'r') as filee:  # edit the mafft file to remove any character which may not work with the perl script
    filedata = filee.read()

# replace spaces to _
remove_characters = [" ", "|", "\\", "/", "^", "(", ")"]
for character in remove_characters:
    filedata = filedata.replace(character, "_")  # replace everyhting in the id names which is a different character
    # with an underscore

# write file out again
with open(file, 'w') as filee:  # save the file as an outfile
    filee.write(filedata)

# READ IN FILE
record_dict = SeqIO.index(file, "fasta")  # read in the AA mafft alignment fasta file

listed = list(record_dict)  # list of ids/records from the dictionary of sequences from the fasta

ids = []
good = 0
new_list = []
new_x = []
for x in listed:
    week = 0  # week equal to zero as no metadata file
    sequence = (record_dict[x].seq)  # the sequence of the associated id/name
    if (startweek <= week) & (week <= endweek):   # all sequences included
        new_x.append(x)  # add all records
        ids.append(x)
        good += 1
        position = 0
        new_list.append(list(sequence))  # create a list of the sequence

print(f"\n {good} sequences from {file} fell between weeks {startweek} and {endweek}\n ")


sequence_table = pandas.DataFrame(index=new_x, data=new_list)  # Create dataframe of the sequences.
sequence_table.columns += 1  # add 1 to the column names as starts from zero, need to start at base 1 = 1

types = len(list(ids))
minii = types * min_maf  # length of array - multiply
if minii < abs_min:   # if minimum is less than the absolute minimum, change to absolute minimum
    minii = abs_min

with open('%s.pipeline.log' % file, 'a') as f: #  print out to log
    f.write(f"\n {good} sequences from {file} fell between weeks {startweek} and {endweek}\n ")  # print to log
    f.write(f"\n Loaded {types} sequences")

with open('%s_covid_variants.tdt' % file, 'w') as g:
    for each_id in ids:
        g.write(f"\t{each_id}")   # for output file, print header with ids used

    g.write(f"\n")

sequence_table.replace(to_replace='[^A-Z\-]', value='Z', regex=True)   # replace any non AA with a Z value
# Replace anything not a letter \ or - with Z.
# For every line/file in the sequence table,

size = sequence_table.shape  # find dimensions of dataframe
final_output = []
snp_count = 0
for i in range(1, size[1]):  # for each column in the length of row 1
    bad = 0
    base_count = []
    col = sequence_table[i]  # column number
    for j in range(0, col.size): # for each row in the column
        if col[j] == None:  # if no entry in the column, continue
            continue
        if re.match('([A-WY])', col[j]):   # if the entry is an AA
            base_count.append(col[j])      # add base to count
            base_count_join = ''.join(base_count)  # remove spaces in list
        else:
            bad = bad + 1  # add to bad count if not an AA
            base_count.append(col[j]) # append to list of
            base_count_join = ''.join(base_count) # remove spaces

    counted = collections.Counter(base_count) # tally the number of counts per base for each column
    new_counted = counted.most_common(3)  # gives the 2 most common AA occurances and their counts in order.
    # therefore the AA in position zero is the Major AA and the value in position one is the minor AA

    if len(new_counted) == 1: # if only one AA counted, not a snp, continue to next column
        continue

    if len(new_counted) == 2:    # if length 2, assign first and second most common to major and minor
        major_base = new_counted[0][0]
        major_count = new_counted[0][1]

        minor_base = new_counted[1][0]
        minor_count = new_counted[1][1]

        if major_base == "-":   # if major - or X, continue
            continue

        if major_base == "X":
            continue

        if minor_base == "-":   # if minor - or X continue
            continue

        if minor_base == "X":
            continue

    if len(new_counted) == 3:      # for those with length 3
        major_base = new_counted[0][0]
        major_count = new_counted[0][1]

        minor_base = new_counted[1][0]
        minor_count = new_counted[1][1]

        if major_base == "-":
            continue

        if major_base == "X":
            continue

        if minor_base == "-":   # for length 3, if minor is - or X, use next most common count
            minor_base = new_counted[2][0]
            minor_count = new_counted[2][1]
            if minor_base == "X":
                continue

        if minor_base == "X":
            minor_base = new_counted[2][0]
            minor_count = new_counted[2][1]
            if minor_base == "-":
                continue

    # if($n > $min && 1-($bad/$types) >$min_call_rate && $banned{$i} <1)
    if (minor_count > minii) and (1 - (bad / types) > min_call_rate):  # if minor base count greater than minimum
        # threshold and enough good counts to be considered a snp
        snp_count += 1  # add one to snp count
        with open('%s_covid_variants.tdt' % file, 'a') as g:   # print snp to file
            g.write("Base " + str(i) + " " + str(major_base) + str(major_count) + str(minor_base) + str(
                minor_count) + "\t" + str(base_count_join) + "\n")


print(
    f"\n Found {snp_count} SNPs. The minor allele frequency used to discriminate the SNPs was {minii}.")
