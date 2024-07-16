from Bio import SeqIO, AlignIO, Align
from Bio.PDB import PDBParser
import numpy as np
from Bio.SeqUtils import seq1 as three2one
import matplotlib as plt
 
#! WORK IN PROGRESS 

#TODO: Load contact matrix from JSON
#TODO: Add residue distance in final output matrix
#TODO: Accept Newick trees for final data package
#TODO: Add output as JSON
#TODO: Add extra error checkpoints
#TODO: Add Customizability (thresholds etc.) and API Calls format
#TODO: Add filtration system to load files
#TODO: Add matplotlib to visualize results
#TODO: Annotate code for repo
#TODO: Reformat + Annotate Code
#TODO: Format output matrix for NN use
#TODO: Format output matrix for JBrowse use 

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
            if distance < distance_threshold:
                contact_matrix[i, j] = contact_matrix[j, i] = 1
    
    return contact_matrix, "".join([three2one(atom.get_parent().get_resname()) for atom in ca_atoms])

# Load TreeFam alignment
def load_treefam_alignment(alignment_file):
    return AlignIO.read(alignment_file, "fasta")

def create_pdb_to_alignment_mapping(pdb_sequence, treefam_alignment):
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

    if best_score > 1.5:  # Adjust this threshold as needed
        pdb_to_aln = {}
        aligned_pdb, aligned_treefam = best_alignment.aligned
        pdb_start = aligned_pdb[0][0]
        treefam_start = aligned_treefam[0][0]
        
        for (pdb_slice, treefam_slice) in zip(aligned_pdb, aligned_treefam):
            for i, j in zip(range(pdb_slice[0], pdb_slice[1]), range(treefam_slice[0], treefam_slice[1])):
                pdb_to_aln[i] = j - treefam_start

        return pdb_to_aln, best_record.id
    else:
        raise ValueError(f"No sufficiently matching sequence found. Best match score: {best_score:.2f}")


# Map contact matrix to alignment
def map_contacts_to_alignment(contact_matrix, pdb_to_aln_mapping, alignment_length):
    alignment_contacts = np.zeros((alignment_length, alignment_length))
    for i in range(contact_matrix.shape[0]):
        for j in range(contact_matrix.shape[1]):
            if contact_matrix[i, j] == 1:
                if i in pdb_to_aln_mapping and j in pdb_to_aln_mapping:
                    aln_i = pdb_to_aln_mapping[i]
                    aln_j = pdb_to_aln_mapping[j]
                    alignment_contacts[aln_i, aln_j] = 1
    return alignment_contacts


def visualize_results(contact_matrix, alignment_contacts, pdb_sequence, matching_seq_id):
    pass

# Main workflow
pdb_file = "data/1jnx.pdb"
treefam_alignment_file = "data/TF105060.fa"
contact_matrix, pdb_sequence = load_pdb_and_generate_contacts(pdb_file)
treefam_alignment = load_treefam_alignment(treefam_alignment_file)
pdb_to_aln_mapping, matching_seq_id = create_pdb_to_alignment_mapping(pdb_sequence, treefam_alignment)
alignment_contacts = map_contacts_to_alignment(contact_matrix, pdb_to_aln_mapping, len(treefam_alignment[0]))

visualize_results(contact_matrix, alignment_contacts, pdb_sequence, matching_seq_id)