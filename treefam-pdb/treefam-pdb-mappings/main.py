#Krishna Bhatt @ Holmes Lab (UC Berkeley) 2024
#Integrated Pipeline Product

#! CONSTRUCTION IN PROGRESS

#*
#* Imports =======================================================================================================
#*

import os
import numpy as np

from Bio import AlignIO, Align
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1 as three2one

from utils import numpy2json, json2numpy

#*
#* Data Loading ==================================================================================================
#*

# Loading TreeFam Alignments with AlignIO Module
def load_treefam_alignment(alignment_file):
    return AlignIO.read(alignment_file, "fasta")

# Loading PDB Files (headers, structures, residue sequences)
def load_pdb_family(fam_directory):
    parser = PDBParser()
    pdb_family = []
    
    for filename in os.listdir(fam_directory):
        pdb_file = os.path.join(fam_directory, filename)
        
        structure = parser.get_structure("protein", pdb_file)
        ca_atoms = [atom for atom in structure.get_atoms() if atom.name == "CA"]
        res_sequence = "".join([three2one(atom.get_parent().get_resname()) for atom in ca_atoms])
        
        pdb_family.append((parser.get_header("protein", pdb_file), structure, res_sequence))

    return pdb_family

# Load Distance Matrices from JSON Files
def load_distances(dist_directory):
    dist_family = []
    
    for filename in os.listdir(dist_directory):
        dist_file = os.path.join(dist_directory, filename)
        distance_matrix = json2numpy(dist_file)
        dist_family.append(distance_matrix)
    
    return dist_family

#*
#* Distance Matrices to Alignment Mapping ========================================================================
#*
       
# Create a dictionary mapping from each PDB residue to corresponding TreeFam columns
def pdb2tree_mapping(pdb_seq, tree_seq, match_threshold=1.5):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    pdbtree_aln = aligner.align(pdb_seq, str(tree_seq.seq))[0]
    score = pdbtree_aln.score / len(pdb_seq)
    
    if score > match_threshold:
        pdb2tree = {}
        aligned_pdb, aligned_treefam = pdbtree_aln.aligned
        tree_start = aligned_treefam[0][0]
        
        for (pdb_slice, tree_slice) in zip(aligned_pdb, aligned_treefam):
            for i, j in zip(range(pdb_slice[0], pdb_slice[1]), range(tree_slice[0], tree_slice[1])):
                pdb2tree[i] = j - tree_start
        
        return pdb2tree, score
    else:
        raise ValueError("PDB Structure Sequence does not sufficiently match the TreeFam Residue Sequence")

# Mapping Distance Matrix to TreeFam Alignment Shape
def dist2tree_mapping(dist_matrix, pdb2tree_map, treeSeq_length):
    tree_contacts = np.full((treeSeq_length, treeSeq_length), -1)
    
    for i in range(dist_matrix.shape[0]):
        for j in range(dist_matrix.shape[1]):
            if i in pdb2tree_map and j in pdb2tree_map:
                aln_i = pdb2tree_map[i]
                aln_j = pdb2tree_map[j]
                tree_contacts[aln_i, aln_j] = dist_matrix[i, j]
                
    return tree_contacts

#*
#* Key Feature Determination =====================================================================================
#*



#*
#* Pipeline Integration & Testing ================================================================================
#*

#TODO: Construct main pipeline based on others code and main runtime file
