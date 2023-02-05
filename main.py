#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 11:39:57 2022
Description of the program available at :
https://github.com/mnicolleUTC/Phylogenetic-Analyses.git
@author: Nicolle Mathieu
Based on https://taylor-lindsay.github.io/phylogenetics/
"""

import os
import Bio as Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator 
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt

def prepare_aligment(file):
     """
     Generate a scala file in which all sequence has been complete by '-' in 
     order to equalize length of all sequence
     If not applied generate ValueError: Sequences must all be the same length
     when using AlignIO.read()

     Args:
         file (str): Name of the input file with different sequences length

     Returns:
         output_file (str): Name of the output file with equalized sequences
         length
     """
     records = SeqIO.parse(input_file, 'fasta')
     records = list(records) # make a copy, otherwise our generator
                             # is exhausted after calculating maxlen
     maxlen = max(len(record.seq) for record in records)
     # pad sequences so that they all have the same length
     for record in records:
          if len(record.seq) != maxlen:
               sequence = str(record.seq).ljust(maxlen, '-')
               record.seq = Seq.Seq(sequence)
     assert all(len(record.seq) == maxlen for record in records)
     # write to temporary file and do alignment
     output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
     with open(output_file, 'w') as f:
          SeqIO.write(records, f, 'fasta')
     return output_file
     
def inspect_file(file):
     """
     Inspect the dataset and print main information of each record

     Args:
         file (str): Input file that contains sequences data
     """
     # Infer sequence from file 
     seq = SeqIO.parse(file, "fasta")
     # Convert to list in order to see content
     list_seq = list(seq)
     for i,elt in enumerate(list_seq):
          print(list_seq[i])
          print(f"Number of element {len(list(list_seq[i]))}\n")
     

if __name__=="__main__":
     # Define path to input data
     input_file = "data.fasta"
     # Inspect file
     #inspect_file(input_file)
     # Prepare for aligment (Sequence must be of same length)
     patted_file = prepare_aligment(input_file)
     # Make alignment
     align = AlignIO.read(patted_file, "fasta")
     # Export alignment to a output_file with phylip format
     AlignIO.write(align, "save_align.phy", "phylip")
     # Read alignment from previous file
     with open("save_align.phy","r") as aln: 
	     alignment = AlignIO.read(aln,"phylip")
     # Instantiate calculator 
     # For this analysis, we will use the ‘identity’ model
     calculator = DistanceCalculator('identity')
     # Generate distance matrix data
     distance_matrix = calculator.get_distance(alignment)
     # Save matrix data to an output text file
     with open("matrice_distance_value.txt",'w') as file:
          file.write(str(distance_matrix))
     # Instantiate tree constructor
     constructor = DistanceTreeConstructor(calculator)
     # Build the tree
     tree = constructor.build_tree(alignment)
     # Plot the tree in terminal
     Phylo.draw_ascii(tree)
