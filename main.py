#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 11:39:57 2022
Description of the program available at :
https://github.com/mnicolleUTC/Phylogenetic-Analyses.git
@author: Nicolle Mathieu
Based on https://taylor-lindsay.github.io/phylogenetics/
Run in conda environment in python 3.9 
"""

import os
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator 
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from ete3 import Tree, TreeStyle
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
    records = SeqIO.parse(file, 'fasta')
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
    output_file = '{}_padded.fasta'.format(os.path.splitext(file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    return output_file

def build_tree_for_ete3(tree_biopython):
    """Convert a biopython tree to an ete3 tree

    Args:
        tree_biopython (Bio.Phylo.BaseTree): Biopython tree builded with input
        data sequences

    Returns:
        builded_tree (ete3.coretype.tree.TreeNode): Same tree in format ete3 
    """
    clade = tree_biopython.root
    builded_tree = Tree()
    for child in clade.clades:
        child_tree = build_tree_for_ete3(child)
        builded_tree.add_child(child=child_tree, name=child.name, dist=child.branch_length)
    return builded_tree
     
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
          
def generate_tree(fasta_file):
    """Generate a tree based on sequence data provided in fasta_file

    Args:
        fasta_file (str): Input file that contains sequences data
        
    Returns:
        tree (Bio.Phylo.BaseTree): Tree object from Bio.Phylo package
    """
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
    return tree
     
     
if __name__=="__main__":
    # Define path to input data
    input_file = "data.fasta"
    # Inspect file
    #inspect_file(input_file)
    # Build tree from input_data
    tree = generate_tree(input_file)
    # Plot the tree in terminal
    #Phylo.draw_ascii(tree)
    ete3_tree = build_tree_for_ete3(tree)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.arc_start = -180 # 0 degrees = 3 o'clock
    ts.arc_span = 180
    ts.show_leaf_name = True
    ete3_tree.show(tree_style=ts)