#Krishna Bhatt @ Holmes Lab (UC Berkeley) 2024
#Standalone testing for new features (Outside of pipeline)

import numpy as np
import os

from Bio import AlignIO, Align
from Bio.PDB import PDBParser, Superimposer
from Bio.SeqUtils import seq1 as three2one

from utils import numpy2json, json2numpy, numpy2csv
 
#! CONSTRUCTION ZONE

#TODO: Develop evolseqpair class type including contacts, alignment scores, function to determine evolutionary divergence (newick trees)
# Feature-wise Evol closeness between sequences
'''
Key Feature Determination:
* Inputs (which are features of a matching DNA segment in an alignment) to function that must exceed threshold to be called key
Inputs- 
    Length of segment
    Protein shape variation of segment (MSE kinda thing)
    Number of Protein contacts in area
'''
# Key Residue Sites and Structures (Scoring System) --> To visualization
# Estimated Key Feature Emergence
#TODO: Add pop up description in matplotlib view (per-residue alignment scores, distance, etc.)
#TODO: Add functionality to Plotly protein models with additional evolutionary trends
#TODO: Annotate Code

# Load PDB structure and generate contact matrix
def load_pdb_and_generate_contacts(pdb_file, distance_threshold=8.0):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    # Extract C-alpha atoms
    ca_atoms = [atom for atom in structure.get_atoms() if atom.name == "CA"]
    
    # Generate contact matrix
    n = len(ca_atoms)
    contact_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            distance = ca_atoms[i] - ca_atoms[j]
            contact_matrix[i, j] = contact_matrix[j, i] = distance

    return contact_matrix, "".join([three2one(atom.get_parent().get_resname()) for atom in ca_atoms])

# Load TreeFam alignment
def load_treefam_alignment(alignment_file):
    return AlignIO.read(alignment_file, "fasta")

def create_pdb_to_alignment_mapping(pdb_sequence, treefam_alignment, match_threshold):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    best_score = 0
    best_alignment = None
    best_record = None

    for record in treefam_alignment:
        alignment = aligner.align(pdb_sequence, str(record.seq))[0]
        score = alignment.score / len(pdb_sequence)
        print(f"Sequence: {record.id}, Score: {score:.2f}")
        if score > best_score:
            best_score = score
            best_alignment = alignment
            best_record = record

    print(f"Best score: {best_score:.2f}")
    print(f"PDB sequence length: {len(pdb_sequence)}")
    print(f"Best matching sequence length: {len(str(best_record.seq))}")

    if best_score > match_threshold:  # Adjust this threshold as needed
        pdb_to_aln = {}
        aligned_pdb, aligned_treefam = best_alignment.aligned
        pdb_start = aligned_pdb[0][0]
        treefam_start = aligned_treefam[0][0]
        
        for (pdb_slice, treefam_slice) in zip(aligned_pdb, aligned_treefam):
            for i, j in zip(range(pdb_slice[0], pdb_slice[1]), range(treefam_slice[0], treefam_slice[1])):
                pdb_to_aln[i] = j - treefam_start

        return pdb_to_aln, best_record.id, best_alignment
    else:
        raise ValueError(f"No sufficiently matching sequence found. Best match score: {best_score:.2f}")


# Map contact matrix to alignment
def map_contacts_to_alignment(contact_matrix, pdb_to_aln_mapping, alignment_length):
    
    alignment_contacts = np.full((alignment_length, alignment_length), -1)
    
    for i in range(contact_matrix.shape[0]):
        for j in range(contact_matrix.shape[1]):
            if i in pdb_to_aln_mapping and j in pdb_to_aln_mapping:
                aln_i = pdb_to_aln_mapping[i]
                aln_j = pdb_to_aln_mapping[j]
                alignment_contacts[aln_i, aln_j] = contact_matrix[i, j]
                
    return alignment_contacts

#Main Mapping workflow
def mapping(alignment_f, pdb_f, contact_matrix_f, out, match_threshold=1.5):
    alignment = AlignIO.read(alignment_f, "fasta")
    cm, pdb_sequence = load_pdb_and_generate_contacts(pdb_f)
    contact_matrix = json2numpy(contact_matrix_f)
    
    pdb2aln, matchingID = create_pdb_to_alignment_mapping(pdb_sequence, alignment, match_threshold)
    alignment_contacts = map_contacts_to_alignment(contact_matrix, pdb2aln, len(alignment[0]))
    numpy2json(alignment_contacts, out)
    return(alignment_contacts)

#Testing Workflow
'''
if __name__ == "__main__":
    
    pdb_file = "treefam-pdb/treefam-pdb-mappings/samples/1jnx.pdb"
    treefam_alignment_file = "treefam-pdb/treefam-pdb-mappings/samples/TF105060.fa"
    contact_matrix, pdb_sequence = load_pdb_and_generate_contacts(pdb_file)
    treefam_alignment = load_treefam_alignment(treefam_alignment_file)
    pdb_to_aln_mapping, matching_seq_id = create_pdb_to_alignment_mapping(pdb_sequence, treefam_alignment, 1.5)
    alignment_contacts = map_contacts_to_alignment(contact_matrix, pdb_to_aln_mapping, len(treefam_alignment[0]))

    print(alignment_contacts.shape)
    numpy2json(alignment_contacts, "mapping.json")
    new = json2numpy("mapping.json")
    print(new.shape)
    numpy2json(new, "mapping2.json")
'''

def genContacts(dist_matrix, dist_thresh=8):
    contacts = np.zeros(dist_matrix.shape)
    for i in range(0, dist_matrix.shape[0]):
        for j in range(0, dist_matrix.shape[1]):
            if dist_matrix[i][j] < dist_thresh:
                contacts[i][j] = contacts[j][i] = 1
    return contacts

def keyContactAreas (dist_matrix, dist_thresh=8, seqDist_thresh=5, numContact_thresh=5, strictness=2):
    score = 0
    contactSegs = [[]]
    
    for i in range(dist_matrix.shape[0]):
        resContacts = 0
        
        for j in range(dist_matrix.shape[1]):
            if dist_matrix[i][j] < dist_thresh and abs(i - j) > seqDist_thresh:
                resContacts += 1
        
        if (resContacts > numContact_thresh):
            score = strictness
        else: score -= 1
                
        if score > 0:
            contactSegs[-1].append(i+1)
        elif len(contactSegs[-1]) != 0:
            contactSegs.append([])
    
    return contactSegs


def keyAlnAreas(alignment, numMatches_thresh=0.9, strictness=2):
    score = 0
    keySegs = [[]]
    
    for i in range(len(alignment[0].seq)):
        resMatches = {}
        passed = False
        
        for j in range(len(alignment)):
            if (alignment[j].seq)[i] in resMatches:
                resMatches[(alignment[j].seq)[i]] += 1
            else:
                resMatches[(alignment[j].seq)[i]] = 1
        
        for res in resMatches:
            if resMatches[res] > numMatches_thresh * len(alignment):
                score = strictness
                passed = True
                
        if not passed: score -= 1
                
        if score > 0:
            keySegs[-1].append(i+1)
        elif len(keySegs[-1]) != 0:
            keySegs.append([])
    
    return keySegs
            
def imposeStructure(pdb_structures):
    ref_structure = pdb_structures[0]
    aligner = Superimposer()
    
    for structure in pdb_structures[1:]:
        ref_atoms = [atom for atom in ref_structure.get_atoms() if atom.name == 'CA']
        atoms = [atom for atom in structure.get_atoms() if atom.name == 'CA']
        
        aligner.set_atoms(ref_atoms, atoms)
        aligner.apply(structure.get_atoms())
    
    return pdb_structures

def calcVariance(imposed_structures):
    all_coords = []
    
    for structure in imposed_structures:
        coords = [atom.coord for atom in structure.get_atoms() if atom.name == 'CA']
        all_coords.append(coords)

    all_coords = np.array(all_coords)
    variance = np.var(all_coords, axis=0)
    
    return np.mean(variance, axis=0)

def keyVarAreas(variance, var_thresh=0.5, strictness=2):
    score = 0
    keySegs = [[]]
    
    for i, var in enumerate(variance):
        if var < var_thresh:
            score = strictness
        else:
            score -= 1
            
        if score > 0:
            keySegs[-1].append(i)
        elif len(keySegs[-1]) != 0:
            keySegs.append([])
    
    return keySegs

if __name__ == '__main__':
    pdb_file = "treefam-pdb/treefam-pdb-mappings/samples/1jnx.pdb"
    treefam_alignment_file = "treefam-pdb/treefam-pdb-mappings/samples/TF105060.fa"
    dist_matrix, pdb_sequence = load_pdb_and_generate_contacts(pdb_file)
    treefam_alignment = load_treefam_alignment(treefam_alignment_file)
    #numpy2csv(genContacts(dist_matrix), 'contacts.csv')
    
    #print(contactArea(dist_matrix))
    print(keyAlnAreas(treefam_alignment))