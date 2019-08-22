#!/share/apps/python-2.7.2/bin/python
#
# Either include pyrosetta in your path or run something like 
# source /home/ateufel/Pyrosetta/PyRosetta.monolith.ubuntu.release-80/SetPyRosettaEnvironment.sh
# so that pyrosetta is in the path
# run this program like: python score_mutant_AIT.py tiny_gene
# where tiny_gene is the name of the pdb file you want to use without the ".pdb" part
# note that your pdb files must start a residue 1, this may mean you have to renumber your pbd file 


# Imports
import sys
import argparse
import random, math, os 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover

from toolbox import mutate_residue
from time import time
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_filename', action="store", type=str)
    parser.add_argument('replicate_number', action="store", type=int)

    inputs = parser.parse_args()    
    #takes name of pdb file without the extention
    pdb_file = inputs.pdb_filename
    prot_name = pdb_file.split('/')[-1].split('.')[0]
    #set up timer to figure out how long the code took to run
    t0=time()
    fasta_file = pdb_file.replace('/structures/', '/fastas/').replace('.pdb', '.fasta')
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    assert len(records)==1
    wt_seq = str(records[0].seq)   
  
    # Initialize Rosetta.
    #init(extra_options='-mute basic -mute core')
    init(extra_options='-mute basic -mute core -rebuild_disulf false -detect_disulf false')

    ########################
    # Constants
    ########################
    PACK_RADIUS = 12.0
    #Amino acids
    AAs = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    AAs_choice_dict = {}
    for aa in AAs:
        AAs_choice_dict[aa] = [other_aa for other_aa in AAs if other_aa != aa]
    #Number of mutations to accept
    max_accept_mut = 10*len(wt_seq)
    #max_accept_mut = 2048

    #Population size
    N = 1000
    #Beta (temp term)
    beta = 1
    #Fraction of the WT stability value to shoot for 
    threshold_fraction = 0.5
    ########################
    ########################

    #Prepare data headers
    data = ['Variant,Rosetta Score,"delta-delta-G",Probability,Generation\n']

    #Load a clean pdb file  
    initial_pose = pose_from_pdb(pdb_file)
    if '.clean' in pdb_file:
        pdb_file = ''.join(pdb_file.split('.clean'))

    #Set up ScoreFunction
    sf = get_fa_scorefxn()

    #Set up MoveMap.
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)

    #Pack and minimize initial pose to remove clashes.
    pre_pre_packing_score = sf(initial_pose)

    task = standard_packer_task(initial_pose)
    task.restrict_to_repacking()
    task.or_include_current(True)
    pack_rotamers_mover = RotamerTrialsMover(sf, task)
    pack_rotamers_mover.apply(initial_pose)

    min_mover = MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(sf)
    min_mover.min_type('dfpmin_armijo_nonmonotone')
    min_mover.apply(initial_pose)

    post_pre_packing_score = sf(initial_pose)
    
    #Threshold for selection 
    threshold = post_pre_packing_score * threshold_fraction
    print 'threshold:', threshold

    data.append('WT,' + str(post_pre_packing_score) + ',0.0,0.0,0\n')

    #number of residues to select from
    n_res = initial_pose.total_residue()

    #start evolution
    i=0
    gen=0
    while i < max_accept_mut:

        #update the number of generations that have pased
        gen+=1

        #print 'accepts:', i 

        #pick a place to mutate
        mut_location = random.randint(1, n_res)

        #get the amino acid at that position
        res = initial_pose.residue(mut_location)
        
        #choose the amino acid to mutate to
        #new_mut_key = random.randint(0,len(AAs)-1)
        #proposed_res = AAs[new_mut_key]
        proposed_res = random.choice(AAs_choice_dict[res.name1()])
 
        #make the mutation
        mutant_pose = mutate_residue(initial_pose, mut_location, proposed_res, PACK_RADIUS, sf)

        #score mutant
        variant_score = sf(mutant_pose)

        #get the probability that the mutation will be accepted
        probability = calc_prob_mh(variant_score, post_pre_packing_score, N, beta, threshold)

        #test to see if mutation is accepted
        if random.random() < probability:
    
            #create a name for the mutant if its going to be kept 
            variant_name = res.name1() + str(initial_pose.pdb_info().number(mut_location)) + str(proposed_res)

            #save name and energy change
            data.append(variant_name + "," + str(variant_score) + "," + str(variant_score - post_pre_packing_score) + "," + str(probability) + "," + str(gen) + "\n")

#            if i == (max_accept_mut - 1):
#                final_pdb_name=pdb_file.replace('.pdb', '_thresh={}_Neff={}_beta={}_i={}_nmut={}.pdb'.format(threshold_fraction, N, beta, inputs.replicate_number, i))  
#                mutant_pose.dump_pdb(final_pdb_name)

            #update the wildtype 
            initial_pose = mutant_pose
            post_pre_packing_score = variant_score

            #update number of accepts
            i+=1

    print '\nMutations and scoring complete.'
    t1 = time()
    # Output results.
    output_filename = '../Results/{}/{}_thresh={}_Neff={}_beta={}_i={}.csv'.format(prot_name, prot_name, threshold_fraction, N, beta, inputs.replicate_number)
    with open(output_filename, "w") as outfile:
        outfile.writelines(data)

    print 'Data written to:', output_filename
    print 'program takes %f' %(t1-t0)


###assorted functions that have to do with scoring and prob of acceptance ####


#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x(stab_org, beta, thresholds)
  xj = calc_x(stab_mut, beta, thresholds)

  if xj > xi:
    return((1.0))
  else:
    exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))



#score functions for met-hastings selection
def calc_prob_fix(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x_fix(stab_org, beta, thresholds)
  xj = calc_x_fix(stab_mut, beta, thresholds)

  if xi == xj:
    return(1/float(N))

  if(xj==0.0):
    return(0.0)

  try:
    p =((1-pow( (xi/xj),2)) /(1-pow( (xi/xj), (2 * float(N) )) ) )
  except OverflowError as e:
    p = 0.0
  return (p)



def calc_x_fix(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total = 1/(safe_calc(exponent) + 1)
  return(total)

def calc_x(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total += -math.log(safe_calc(exponent) + 1)
  return(total)


def safe_calc(exponent):
  if exponent > 700:
    #print("system maxed")
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))


#Run main program
if __name__ == '__main__':
    main()
#    test_function()
