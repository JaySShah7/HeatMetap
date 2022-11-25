#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 11:24:51 2022

@author: jay
"""
import subprocess
import pandas as pd
import os
import gzip
import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import folium
from folium.plugins import HeatMap

def count_residues(filename):       
    # parse .gz files
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            fasta=f.readlines()
    # else assume uncompressed fasta file
    else:
        with open(filename,'r') as f:
            fasta=f.readlines()
                
    
    # number of sequences is number of headers
    # number of header is number of line starting with '>'
    num_sequences = len([x for x in fasta if x.startswith('>')])
    if not num_sequences == 1:
        raise Exception("More than 1 sequence in reference file!!!")
        
    
    # Get all non-header i.e. residue lines
    residues = [x for x in fasta if not x.startswith('>')]  
    num_residues = 0
    
    for residue in residues:
        # add length of residue line, removing the newline character
        # NOTE: We should also remove the '*' i.e. stop codon? as its not a residue?
        num_residues += len(residue.replace(' ', '').replace('\n', ''))
        
    return num_residues


# Function which parses list of MDS strings and outputs list of percent identity
def parse_md(mds):
    digits = "(\d+)"
    nondigits = "(\D+)"
    percent_identities = []
    
    
    for md in mds: 
        matches  = 0
        mismatches = 0
        # Sum up all matches
        matches_list = re.findall(digits, md)
        for match in matches_list:
            matches += int(match)
            
        # Sum up all mismatches
        mismatches_list = re.findall(nondigits, md)
        for mismatch in mismatches_list[1:]:
            mismatches+=len(mismatch)

        # calculate percent identity
        identity = 100 * matches / (matches+mismatches)
        percent_identities.append(identity)
    return percent_identities


# function which takes SAM filename as input and outputs 3 lists: read positions,
# CIGAR strings, and MD tags
def parse_sam(SAM):
    cigars = []
    mds = []
    positions = []
    lengths = []
    # Parse SAM to get cigar strings
    with open(SAM, 'r') as f:
        for line in f: 
            # Lines starting with @ are headers
            if not line.startswith("@"):
                # cigar string is 5th element of line
                line_columns = line.rstrip().split()
                # CIGAR string is column 5
                CIGAR = line_columns[5]
                # In case unmapped read
                if CIGAR == "0":
                    continue
                cigars.append(CIGAR)
                # Detect MD string - MD string is optional
                try:
                    md = [x for x in line_columns if x.startswith("MD")][0]
                    mds.append(md)
                except IndexError:
                    pass
                
                # Extract position
                positions.append(int(line_columns[3]))
                
                # Extra data
                seq_len = (len(line_columns[9]))
                lengths.append(seq_len)
                
                
        return positions, cigars, mds, lengths

# Function which takes list of read positions, read lengths, read percent identities
# as input and outputs overall percent genome coverage in sample as well as median
# coverage per base
# This method is not 100% accurate, but by limiting to percent identity >50% it 
# provides a good estimate of genome coverage
def calculate_coverage(positions, lengths, percent_identities, min_coverage=1,
                       min_percent_identity=50, reference_length=0):
    # Create array of zeros representing each bases coverage
    if reference_length:
        reference = np.zeros(reference_length)
    else: # if reference length not given, calculate best approximation
        reference = np.zeros(max(positions) + lengths[-1])    
    
    # go through each base
    for i, position in enumerate(positions):
        # only count as mapped if PI > 50%
        if percent_identities[i] >= min_percent_identity:
            # add 1 to each base covered
            for j in range(lengths[i]):
                if position+j-1 >=reference_length:
                    continue
                reference[position+j-1]+=1
    # return count of bases above threshold / total bases i.e. percentage
    # and also median coverage i.e. median counts
    return 100* np.count_nonzero(reference>=min_coverage)/len(reference), np.median(reference)


#%%
if __name__ == "__main__":
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    # Required arguments
    parser.add_argument("input", help="input folder containing read data (input_data.csv)")
    parser.add_argument("output", help="desired output folder")
    parser.add_argument("reference", help="reference genome to track")
    
    # Optional arguments
    parser.add_argument("-d", "--download", help="[EXPERIMENTAL] attempt to automatically download read files into input folder")
    parser.add_argument("--cpu", help="number of CPU threads to use (default 12)")
    parser.add_argument("-c", "--coverage", help="minimum coverage of the reference genome for a sample to be included (default 80)")
    parser.add_argument("-i", "--identity", help="minimum percent identity for a mapped read to be included (default 50)")
    parser.add_argument("-r", "--radius", help="output heatmap radius of datapoint (default 60)")
    parser.add_argument("-b", "--blur", help="outmput heatmap blur value per datapoint (default 40)")
    
    # Store arguments to variables
    args = parser.parse_args()
    input_folder = args.input
    output_folder = args.output
    reference = args.reference
    
    # FOR TESTING
    # reference = "Input/Reference/1283077.fasta"
    # input_folder = "Input"
    # output_folder = "Output"
    # cpu=12
    # minimum_identity = 50
    # radius = 70
    # blur = 40
    
    
    if args.cpu:
        cpu = int(args.cpu)
    else:
        cpu = 12
    if args.identity:
        minimum_identity = int(args.identity)
    else:
        minimum_identity  = 50
        if args.coverage:
            minimum_coverage = int(args.coverage)
        else:
            minimum_coverage  = 80  
    if args.radius:
        radius = int(args.radius)
    else:
        radius = 40

    if args.blur:
        blur = int(args.blur)
    else:
        blur = 70
            
    
    #os.chdir(r"/home/jay/JHU/Metagenomics/Project/")
    input_csv = os.path.join(input_folder, "input_data.csv")
    # reference = "Input/Reference/1283077.fasta"
    # cpu = 12
    
    # count length of reference sequence for future
    reference_length = count_residues(reference)
    
    data = pd.read_csv(input_csv)
    sra_numbers = list(data['SRR'])
    
    
    if args.download:
        #this will download the .sra files to output folder
        for sra_id in sra_numbers:
            print ("Currently downloading: " + sra_id)
            prefetch = "prefetch " + sra_id
            print ("The command used was: " + prefetch)
            subprocess.call(prefetch, shell=True)
        
        # this will extract the .sra files from above into a folder named 'Input'
        for sra_id in sra_numbers:
            print ("Generating fastq for: " + sra_id)
            fastq_dump = "fastq-dump --outdir {} --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {}.sra".format(input_folder, sra_id)
            print ("The command used was: " + fastq_dump)
            subprocess.call(fastq_dump, shell=True)
        
        
    #%%
    
    # Build indexes and align 
    
    
    #bowtie2-build Input/Reference/1283077.fasta Input/Reference/1283077
    
    # Create name for index by removing reference file extension, and index using 
    # bowtie2-build
    regex = r"(\S+)\.\S+"
    index_name = re.findall(regex, reference)[0]
    
    # Only file name without extension or directory names (for saving output later)
    index_filename = os.path.basename(index_name)
    subprocess.call("bowtie2-build {} {}".format(reference, index_name), shell=True)
    
    # Now align for each sample
    
    # list to store sam file locations
    output_sams = []
    
    # dictionary to use for heatmap
    heatmap_metrics = {}
    
    for sra in sra_numbers: 
        # find paired reads for sra
        regex = sra + r"\S+"
        # search all files in the Input directory via a regex search
        read_files = re.findall(regex, "   ".join(os.listdir(input_folder)))
        read_files.sort()
        # add input folder name
        read_files = [os.path.join(input_folder, f) for f in read_files]
        
        # output SAM file name
        SAM =  os.path.join(output_folder, "{}_{}.sam".format(index_filename, sra))
        output_sams.append(SAM)
        
        # The total command

        command = r"bowtie2 --no-unal -p {} --very-fast-local -x {} -1 {} -2 {} -S {}".format(
            cpu, index_name, read_files[0], read_files[1], SAM)
        
        print(command)
        if not (os.path.exists(SAM)):
            subprocess.call(command, shell=True)
        
        # Parse SAM file
        positions, cigars, mds, lengths = parse_sam(SAM)
        
        # If there are mapped reads
        if mds:
            percent_identities = parse_md(mds)
        else:
            percent_identities = [0] * len(mds)
        
        # Calculate % coverage of genome, and median coverage per base
        coverage, median_coverage = calculate_coverage(positions, lengths,
                                                       percent_identities,
                                                       reference_length=reference_length,
                                                       min_percent_identity=minimum_identity)
                
        # FRP plot
        plt.scatter(positions, percent_identities, s=1, c=np.array(percent_identities), cmap='rainbow')
        plt.suptitle("Bowtie2 FRP {}".format(SAM))
        plt.title(" Genome Coverage: {:.1f}   Median Coverage: {:.1f}".format(coverage, median_coverage))
        plt.xlabel("Genome Position")
        plt.ylabel("Percent Identity")
        plt.show()
        
        # count number of basepairs in read file for abundance metric
        # in case basepairs present in input CSV
        if 'basepairs_count' in data.columns:
            total_bp = data[data['SRR']==sra]['basepairs_count'].values[0]
        # else count manually
        else:
            total_bp = count_residues(read_files[0])
        # save metric to heatmap metric dictionary
        if coverage < minimum_coverage:
            coverage = 0
        heatmap_metrics[sra] = median_coverage* coverage / total_bp
        
        
    #%% Creation of heatmap

        
    # Scale abundance metric for ideal heatmap
    while min([ x for x in heatmap_metrics.values() if x!=0] ) < 1:
        for sra in heatmap_metrics:
            heatmap_metrics[sra] = heatmap_metrics[sra]*10
    
    heatmap_rows = []
    for sra in sra_numbers:
        # get latitude, longitude
        latitude, longitude = data[data['SRR']==sra][['latitude', 'longitude']].values[0]
        score = heatmap_metrics[sra]
        
        heatmap_row = (latitude, longitude, score)
        heatmap_rows.append(heatmap_row)
    
        
        
    # Calculate center of heatmap using average of all points
    map_center = np.mean([r[0] for r in heatmap_rows]), np.mean([r[1] for r in heatmap_rows])
    
    # Create new heatmap
    hmap = folium.Map(location=map_center, zoom_start=19)
    

    hm= HeatMap(heatmap_rows, min_opacity=0.2, radius=radius, blur=blur,
                max_zoom=1)
    hmap.add_child(hm)
    
    # Add markers
    for row in heatmap_rows:
        folium.Marker(location=[row[0], row[1]]).add_to(hmap)    
    
    hmap.save(os.path.join(output_folder, 'heatmap.html'))


    
    
