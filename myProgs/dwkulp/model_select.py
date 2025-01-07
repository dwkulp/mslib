import os
import subprocess
import multiprocessing
import time
from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
import pandas as pd
from plotnine import ggplot, aes, geom_point, labs

# SELECT: A31_E 3.19,V 18.41 12641,L 15.01 10305,I 13.35 9166
def parse_output_freq_select(file_path):
    mutations = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("SELECT:"):
                parts = line.strip().split(",")
                base_mut = parts[0].strip().split()[1].replace('_', '')[:-1]  # A31_E -> A31
                for part in parts[1:]:
                    mutation = part.split()[0] # V
                    mutations.append(base_mut + mutation)
    return mutations


def find_neighbors(pdb_file, mutation):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    chain_id = mutation[0]
    resnum = int(mutation[1:-1])
    atom_list = []

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[1] == resnum:
                        if 'CB' in residue:
                            atom = residue['CB']
                        elif 'CA' in residue:
                            atom = residue['CA']
                        else:
                            return []
                        atom_list.append(atom)
                        break
                    
    for model in structure:
        for chain in model:
            for residue in chain:
                if chain.id == chain_id and residue.id[1] == resnum:
                    continue
                if 'CB' in residue:
                    atom = residue['CB']
                elif 'CA' in residue:
                    atom = residue['CA']
                else:
                    continue
                atom_list.append(atom)
                
                
    #for atom in atom_list:
    #   print(f"Atom: {atom.get_name()}, Residue number: {atom.get_parent().id[1]}")
    neighbor_search = NeighborSearch(atom_list)
    center_atom = atom_list[0]
    neighbors = neighbor_search.search(center_atom.coord, 6.0)  # 6.0 Å radius for neighbors
    
    neighbors_sorted = sorted(neighbors, key=lambda atom: np.linalg.norm(atom.coord - center_atom.coord))
    
    neighbor_residues = []
    for neighbor in neighbors_sorted:
        #print(f"Neighbor: {neighbor.get_name()}, Residue number: {neighbor.get_parent().id[1]}")
        if neighbor.get_parent().id[1] != resnum:  # Exclude the mutation residue itself
            neighbor_residues.append(neighbor.get_parent())

    return neighbor_residues[:3]  # Return closest 3 neighbors

def create_resfile(mutation, resfile_path, pdb_file):
    neighbors = find_neighbors(pdb_file, mutation)
    #print(f"Neighbors for {mutation}: {neighbors}")
    with open(resfile_path, 'w') as resfile:
        resfile.write("NATRO\nstart\n")
        chain = mutation[0]
        resnum = mutation[1:-1]
        resfile.write(f"{resnum} {chain} PIKAA {mutation[-1]} EX 1 EX 2 EX_CUTOFF 4\n")
        
        for neighbor in neighbors:
            neighbor_chain = neighbor.get_parent().id
            neighbor_resnum = neighbor.id[1]
            resfile.write(f"{neighbor_resnum} {neighbor_chain} NATAA EX 1 EX_CUTOFF 4\n")
            
    resfile_path_des = resfile_path.replace('.txt', '_des.txt')
    with open(resfile_path_des, 'w') as resfile_des:
        resfile_des.write("NATRO\nstart\n")
        resfile_des.write(f"{resnum} {chain} PIKAA {mutation[-1]} EX 1 EX 2 EX_CUTOFF 4\n")
        
        for neighbor in neighbors:
            neighbor_chain = neighbor.get_parent().id
            neighbor_resnum = neighbor.id[1]
            resfile_des.write(f"{neighbor_resnum} {neighbor_chain} ALLAA EX 1 EX 2 EX_CUTOFF 4\n")
    
    resfile_path_wt = resfile_path.replace('.txt', '_wt.txt')
    with open(resfile_path_wt, 'w') as resfile:
        resfile.write("NATRO\nstart\n")
        chain = mutation[0]
        resnum = mutation[1:-1]
        resfile.write(f"{resnum} {chain} NATAA EX 1 EX_CUTOFF 4\n")
    
        for neighbor in neighbors:
            neighbor_chain = neighbor.get_parent().id
            neighbor_resnum = neighbor.id[1]
            resfile.write(f"{neighbor_resnum} {neighbor_chain} NATAA EX 1 EX_CUTOFF 4\n")
    
        

def run_rosetta_docker_fixbb(pdb_file, resfile_path, mutation):
    output_pdb_file = f"{mutation}_fixbb_"
    command = (f"docker run --rm -v {os.getcwd()}:/data rosettacommons/rosetta fixbb -s /data/{pdb_file} -resfile /data/{resfile_path} -nstruct 1 -out:prefix /data/fixbb/{output_pdb_file}")
    print(f"Running command: {command}")
    subprocess.run(command, shell=True)
    full_output_pdb_file = "fixbb/"+output_pdb_file + pdb_file.replace('.pdb', '') + '_0001.pdb'
    return full_output_pdb_file

def run_rosetta_docker_relax(pdb_file, resfile_path, mutation):
    output_pdb_file = f"{mutation}_relax_"
    command = (f"docker run --rm -v {os.getcwd()}:/data rosettacommons/rosetta relax -s /data/{pdb_file} -resfile /data/{resfile_path} -nstruct 1 -relax:quick  -relax:constrain_relax_to_start_coords  -relax:coord_constrain_sidechains -out:prefix /data/relax/{output_pdb_file}")
    print(f"Running command: {command}")
    subprocess.run(command, shell=True)
    full_output_pdb_file = "relax/"+output_pdb_file + pdb_file.replace('.pdb', '') + '_0001.pdb'
    return full_output_pdb_file
    
def run_rosetta_docker_score(pdb_file, mutation):
    command =  (f"docker run --rm -v {os.getcwd()}:/data rosettacommons/rosetta score -s /data/{pdb_file} -out:file:score_only  -out:file:scorefile /data/scores/{mutation}.sc ")
 
    print(f"Running command: {command}")
    subprocess.run(command,shell=True)

def worker(mutation, pdb_file):
    single_resfile_path_fix = f"resfiles/resfile_{mutation}.txt"
    single_resfile_path_des = f"resfiles/resfile_{mutation}_des.txt"
    single_resfile_path_wt  = f"resfiles/resfile_{mutation}_wt.txt"
    create_resfile(mutation, single_resfile_path_fix, pdb_file)
    
    designed_pdb_file = run_rosetta_docker_fixbb(pdb_file, single_resfile_path_fix, mutation)
    designed_pdb_file_relax = run_rosetta_docker_relax(designed_pdb_file, single_resfile_path_fix, mutation)
    run_rosetta_docker_score(designed_pdb_file_relax, mutation)
    
    designed_pdb_file = run_rosetta_docker_fixbb(pdb_file, single_resfile_path_des, mutation+"_des")
    designed_pdb_file_relax = run_rosetta_docker_relax(designed_pdb_file, single_resfile_path_des, mutation+"_des")
    run_rosetta_docker_score(designed_pdb_file_relax, mutation+"_des")
    
    wt_pdb_file = run_rosetta_docker_fixbb(pdb_file, single_resfile_path_wt, mutation+"_wt")
    wt_pdb_file_relax = run_rosetta_docker_relax(wt_pdb_file, single_resfile_path_wt, mutation+"_wt")
    run_rosetta_docker_score(wt_pdb_file_relax, mutation+"_wt")
    
def createMutationScoreFile(mutations, output_freq_select_file):
    for mutation in mutations:
        freq_dict = {}
        with open(output_freq_select_file, 'r') as file:
            for line in file:
                if line.startswith("SELECT:"):
                    parts = line.strip().split(",")
                    base_mut = parts[0].strip().split()[1].replace('_', '')[:-1]  # A31_E -> A31
                    for part in parts[1:]:
                        mutation = part.split()[0]  # V
                        freq = float(part.split()[1])  # Frequency
                        full_mutation = base_mut + mutation
                        freq_dict[full_mutation] = freq
        score_file_wt = f"scores/{mutation}_wt.sc"
        score_file_des = f"scores/{mutation}_des.sc"
        score_file = f"scores/{mutation}.sc"

        if not os.path.exists(score_file_wt) or not os.path.exists(score_file_des) or not os.path.exists(score_file):
            continue

        with open(score_file, 'r') as sc_file, open(score_file_wt, 'r') as wt_file, open(score_file_des, 'r') as des_file:
            wt_score = None
            des_score = None
            aa_score = None

            for line in sc_file:
                if line.startswith("SCORE:"):
                    aa_score = float(line.split()[1])  # Assuming total_score is the second column

            for line in wt_file:
                if line.startswith("SCORE:"):
                    wt_score = float(line.split()[1])  # Assuming total_score is the second column

            for line in des_file:
                if line.startswith("SCORE:"):
                    des_score = float(line.split()[1])  # Assuming total_score is the second column

            if wt_score is not None and des_score is not None and aa_score is not None:
                score_diff1 = des_score - wt_score
                score_diff2 = aa_score - wt_score
                freq = freq_dict.get(mutation, 0)

                with open("mutation_scores.txt", 'a') as output_file:
                    output_file.write(f"{mutation}, {score_diff1}, {score_diff2}, {freq}\n")
                        
def main():
    
    # Make directory structure for this code
    directories = ["resfiles", "scores", "fixbb", "relax"]
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    # Inputs (should be from command line?)
    pdb_file = "5K6I_RSV_AA.pdb"  # Replace with your PDB file name
    output_freq_select_file = "output_freq_select.dat"
    #output_freq_select_file = "test"
    
    #positions_to_include = [26,28,30,32,33,36,37,38,39,40,41,44,46,47,48,49,55,56,57,58,59,61,76,78,79,82,83,86,89,93,96,101,102,103,143,144,145,146,147,148,151,154,158,164,167,171,179,187,188,189,195,198,199,203,207,215,220,221,223,227,230,231,232,233,234,235,236,237,238,239,240,241,242,244,245,246,247,251,252,257,260,261,264,274,277,278,280,281,282,283,284,285,287,288,289,290,296,297,298,299,300,301,302,303,304,306,313,316,318,319,320,321,322,331,332,333,334,336,337,338,339,340,341,342,343,345,349,350,351,352,353,354,355,356,357,358,359,365,366,367,368,369,370,371,372,373,374,375,382,386,391,393,395,397,400,408,410,411,412,413,414,415,416,417,422,423,424,435,441,444,447,450,452,455,458,459,474,475,483,484,485,486,487,488,489,490,491,492,494,495]  # List of positions to include
    positions_to_include = [32,33,47,49]
    
    all_mutations = parse_output_freq_select(output_freq_select_file)
    if len(positions_to_include) == 0:
        mutations = all_mutations
    else:
        mutations = [mutation for mutation in all_mutations if int(mutation[1:-1]) in positions_to_include]
        
    #for mutation in mutations:
    #    print(mutation)
    #    single_resfile_path = f"resfiles/resfile_{mutation}.txt"
    #    create_resfile(mutation, single_resfile_path, pdb_file)
        
    #exit(0)
    
    num_cores = int(multiprocessing.cpu_count() * 0.75)
   

    for i in range(0, len(mutations), num_cores):
        batch = mutations[i:i + num_cores]
        pool = multiprocessing.Pool(processes=num_cores)
        pool.starmap(worker, [(mutation, pdb_file) for mutation in batch])
        pool.close()
        pool.join()
        pool = multiprocessing.Pool(processes=num_cores)
        
    # wait for all processes to finish
    print("Waiting for all processes to finish")
    pool.close()
    print("All processes finished")

    # Create a file with mutation scores and frequencies
    createMutationScoreFile(mutations, output_freq_select_file)
    
    # Load the mutation scores into a DataFrame
    mutation_scores = pd.read_csv("mutation_scores.txt", header=None, names=["mutation", "score_diff1", "score_diff2", "freq"])

    # Plot 1: freq vs score_diff1
    plot1 = (ggplot(mutation_scores, aes(x='freq', y='score_diff1')) +
             geom_point() +
             labs(title='Frequency vs Score Difference 1', x='Frequency', y='Score Difference 1'))
    plot1.save("freq_vs_score_diff1.png")

    # Plot 2: freq vs score_diff2
    plot2 = (ggplot(mutation_scores, aes(x='freq', y='score_diff2')) +
             geom_point() +
             labs(title='Frequency vs Score Difference 2', x='Frequency', y='Score Difference 2'))
    plot2.save("freq_vs_score_diff2.png")
        
    print("Done")
        
if __name__ == "__main__":
    main()