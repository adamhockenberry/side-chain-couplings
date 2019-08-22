#import sys
import argparse
#import random, math, os 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover
from Bio.PDB import PDBParser, PDBIO, Dice
import glob
#from toolbox import mutate_residue
#from time import time
#from Bio import SeqIO

init(extra_options='-mute basic -mute core -ignore_zero_occupancy false -rebuild_disulf false -detect_disulf false')

for pdb_file in glob.glob('../Data/structures/*.pdb'):
    print(pdb_file)
    if '.rosetta' in pdb_file:
        continue    
    pdb_file_clean = pdb_file.replace('.pdb', '.rosetta.pdb')
        

    initial_pose = pose_from_pdb(pdb_file)
    initial_pose.dump_pdb(pdb_file_clean)

    io = PDBIO()
    pdb = PDBParser().get_structure(pdb_file.split('/')[-1].split('.')[0], pdb_file_clean)
    chains = list(pdb.get_chains())
    assert len(chains) == 1
    residues = list(chains[0].get_residues())
    Dice.extract(pdb, chains[0].get_id(), 1, len(residues)+1, pdb_file_clean)
    #io.set_structure(chains[0])
    #io.save(pdb_file_clean)
    #io = PDBIO()                

    ##Set up ScoreFunction
    #sf = get_fa_scorefxn()

    ##Set up MoveMap.
    #mm = MoveMap()
    #mm.set_bb(True)
    #mm.set_chi(True)

    ##Pack and minimize initial pose to remove clashes.
    #pre_pre_packing_score = sf(initial_pose)

    #task = standard_packer_task(initial_pose)
    #task.restrict_to_repacking()
    #task.or_include_current(True)
    #pack_rotamers_mover = RotamerTrialsMover(sf, task)
    #pack_rotamers_mover.apply(initial_pose)

    #min_mover = MinMover()
    #min_mover.movemap(mm)
    #min_mover.score_function(sf)
    #min_mover.min_type('dfpmin_armijo_nonmonotone')
    #min_mover.apply(initial_pose)

    #post_pre_packing_score = sf(initial_pose)

    #if i == (max_accept_mut - 1):
    #    final_pdb_name=pdb_file.replace('.pdb', '_thresh={}_Neff={}_beta={}_i={}_nmut={}.pdb'.format(threshold_fraction, N, beta, inputs.replicate_number, i))  
    #    mutant_pose.dump_pdb(final_pdb_name)

 
