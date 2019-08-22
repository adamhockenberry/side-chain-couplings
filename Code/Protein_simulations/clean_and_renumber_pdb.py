#!/share/apps/python-2.7.2/bin/python


# Imports
import argparse
from toolbox import cleanATOM
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa, three_to_one
import glob

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_directory', action="store", type=str)
    inputs = parser.parse_args()    
    #takes name of pdb file without the extention
    for pdb_file in glob.glob(inputs.pdb_directory+'*.pdb'):
        clean_pdb_file = pdb_file.replace('.pdb', '.clean.pdb')
        print('#######################')
        print('#######################{}'.format(pdb_file))
        if 'clean' in pdb_file:
            print('Will overwrite an existing clean pdb so am skipping')
            continue
        
        fasta_outfile_loc = pdb_file.replace('/PDBs/', '/wt_fastas/').replace('.pdb', '.fasta')
        
        #Load and clean up pdb file  
        cleanATOM(pdb_file)
        
        with open(clean_pdb_file, 'r') as infile:
            old_lines = infile.readlines()

        pdb_io = PDB.PDBIO()
        pdb_parser = PDB.PDBParser()
        structure = pdb_parser.get_structure(" ", clean_pdb_file)
        
        if len(structure) != 1:
            print('THERE APPEARS TO BE MORE THAN ONE MODEL IN THIS STRUCTURE BEHAVIOR OF PRORAM IS UNKNOWN ({}). EXITING'.format(clean_pdb_file))
            continue

        chain_counts = {}
        for model in structure:
            for chain in model:
                new_number = 1
                for i, residue in enumerate(chain.get_residues()):
                    res_id = list(residue.id)
                    if res_id[1] != new_number:
                        res_id[1] = new_number
                        residue.id = tuple(res_id)
                    new_number+=1
                chain_counts[chain.id] = new_number
        
        chains = sorted(chain_counts.items(), key=lambda x: x[1])
        chain_to_keep = chains[-1][0]
        chains_to_delete = chains[:-1]
        chains_to_delete = [i for i,j in chains_to_delete]    
        for i,j in enumerate(chains_to_delete):
            structure[0].detach_child(chains_to_delete[i])
        pdb_io.set_structure(structure)
        pdb_io.save(clean_pdb_file)

        
        for model in structure:
            for chain in model:
                print('kept ID {} and deleted {}'.format(chain.id, chains_to_delete))
                seq_list = []
                chainID = chain.get_id()
                for residue in chain:
                    if is_aa(residue.get_resname(), standard=True):
                        seq_list.append(three_to_one(residue.get_resname()))
                    else:
                        seq_list.append('X')
                wt_seq = ''.join(seq_list)

        with open(fasta_outfile_loc, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format('WT', wt_seq))


#Run main program
if __name__ == '__main__':
    main()
