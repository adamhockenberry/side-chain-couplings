# side-chain-couplings
This repository is for the PeerJ manuscript titled [Evolutionary couplings detect side-chain interactions](https://peerj.com/articles/7280/).

The most important thing to note is that the `Data` folder is not included here due to size restrictions. So in order to do more than look at / read through the code, users are required to [download the data archive from Zenodo](https://dx.doi.org/10.5281/zenodo.2552779). 

There might be a few minor modifications to the code that will be necessary just to get directory structures to pair up nicely after you download and unzip the data. My standard format is to have a `Code` (already here), `Data` and `Results` directory here in the project home so all code is written and presumed to follow this basic structure.

The manuscript should resolve any technical questions and this code should be complete, but if users have any issues replicating results they should feel free to contact me via email.

Briefly, for the empirical analyses presented in the paper, users should first turn to the `Code/Empirical_coupling_analyses/` directory. This folder contains a collection of scripts that were used to compile the datasets/alignments using a local install of HMMER and the [representative proteome databases](https://proteininformationresource.org/rps/). Note that the scripts contained in this folder were run on a >100 node cluster with >1TB memory so users should take care if they are repeating any of these analyses. Additionally, the directory names might not line up super nice here since these scripts were run off-line on a separate computer and I didn't take the time to fully edit things into a cohesive directory structure. That said, and without too much detail, homologs were aligned with MAFFT, and the resulting alignemnts were cleaned up and pushed around a bit using these scripts. Then evolutionary coupling analyses were performed separately with `ccmpred_all.sh`, `plmc_all.sh`, and `psicov_all.sh`.

The protein simulation portion was a bit more complicated and was run off line on a local install of ROSETTA, again on an offline cluster so the caveats above apply. Scripts for this portion are contained in `Code/Protein_simulations/`. Of note is that `basic_simulation_include_cys.py` is the main control file here to run single PYROSETTA simulations from input `.PDB` files. The output is a `.csv` of mutations so these are then reconstructed into `.fasta` files and ultimately `.flat` files to be used as input into evolutionary coupling analyses.

Those are the main data generation scripts, which as I noted will be a bit more complicated to run. The remainder of scripts in the `Code` directory should be more straightforward as analysis scripts that work mainly on `.PDB` files (to get contacts) and coupling file outputs from the various programs. These are the notebooks that I feel will be most useful to users and have been annotated and treated with a bit more care for reproducibility and readability. 

But as always, reach out if anything doesn't make sense and/or you need help performing any portion of the data generation / analysis pipeline. At the moment, however, there may not be any more substantial code commits to this repository unless I get feedback looking for clarification.

--Adam
