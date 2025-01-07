import os
import multiprocessing
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys

def get_80_percent_cores():
    total_cores = multiprocessing.cpu_count()
    cores_to_use = max(1, int(total_cores * 0.8))  # Ensure at least 1 core is used
    return cores_to_use

def run_aaScan(pdb, database, skip_existing=True, selection=None, output_file=None):
    # aaScan --pdb ../5K6I_RSV_AA.pdb --fragdb ~/databases/cullpdb/CAonly100/xaa.db 
    # Make output_file that is the name of database file with .dat extension
    skip_it = False
    if output_file is None:
        if selection:
            res_range = selection.split()[-1].replace("resi ", "").replace("-", "_")
            output_file = f"{os.path.splitext(os.path.basename(database))[0]}_aascan_{res_range}.dat"
        else:
            output_file = os.path.splitext(os.path.basename(database))[0] + "_aascan.dat"
    if skip_existing and os.path.exists(output_file):
        skip_it = True
        
    if not skip_it:
        command = f"aaScan --pdb {pdb} --fragdb {database} --sel \"{selection}\" --speak aa_scan --speak PDBFragments --rmsd 0.25 --maxMatches 5000 | grep -E \"DATA|COUNT\" > {output_file}"
        try:
            result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
            print(f"Command output for {database}:", result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error executing command for {database}:", e.stderr)

def main(pdb_file, databases, selection_statements):
    cores = get_80_percent_cores()
    print(f"Running with {cores} concurrent processes...")

    # Use ThreadPoolExecutor to run up to `cores` commands at the same time
    with ThreadPoolExecutor(max_workers=cores) as executor:
        # Start the tasks
        futures = []
        for db in databases:
            for selection in selection_statements:
                output_file = f"{os.path.splitext(os.path.basename(db))[0]}_{selection.replace(' ', '_').replace('and', '').replace('resi', 'res')}.dat"
                futures.append(executor.submit(run_aaScan, pdb_file, db, True, selection, output_file))
        
        # Collect results as they complete
        for future in as_completed(futures):
            try:
                future.result()  # Trigger exception handling if an error occurred
            except Exception as e:
                print("An error occurred:", e)

    print("All tasks completed.")
    
    # list of all output files
    output_files = [f"{os.path.splitext(os.path.basename(db))[0]}_{selection.replace(' ', '_').replace('and', '').replace('resi', 'res')}.dat" for db in databases for selection in selection_statements]
                    
    # For each output file, parse the data, by comma, the first token split by whitespace then the third token is the key, for each other token it is Amino Acid Count Total, sum Count Total for each amino acid
    # An example line is: COUNTS: A91 T   5001	,A 806 5001,L 520 5001,V 383 5001,E 382 5001,K 322 5001,I 278 5001,Q 263 5001,D 249 5001,R 248 5001,S 224 5001,T 210 5001,G 195 5001,F 175 5001,N 167 5001,M 137 5001,Y 135 5001,H 95 5001,W 77 5001,C 74 5001,P 60 5001,X 1 5001
    amino_acid_counts = {}
    total_counts = {}
    for file in output_files:
        with open(file, "r") as f:
            for line in f:
                if line.startswith("COUNTS:"):
                    tokens = line.split(",")
                    key = tokens[0].split()[1] +"_"+ tokens[0].split()[2]
                    for token in tokens[1:]:
                        amino_acid, count, total = token.split()
                        #print(f"key: {key} amino_acid: {amino_acid} count: {count} total: {total}")
                        if key not in amino_acid_counts:
                            amino_acid_counts[key] = {}
                            total_counts[key] = {}
                        amino_acid_counts[key][amino_acid] = amino_acid_counts[key].get(amino_acid, 0) + int(count)
                        total_counts[key][amino_acid] = total_counts[key].get(amino_acid, 0) + int(total)
    
    # Write out the same format as the input files, but with amino acid counts and total counts
    with open("output_count.dat", "w") as f:
        for key in amino_acid_counts:
            f.write(f"COUNTS: {key}")
            for amino_acid in sorted(amino_acid_counts[key], key=lambda aa: amino_acid_counts[key][aa], reverse=True):
                count = amino_acid_counts[key][amino_acid]
                total = total_counts[key][amino_acid]
                f.write(f",{amino_acid} {count} {total}")
            f.write("\n")
            
    # Write out a DATA format with frequencies as the following format:
    # DATA: A   A,28 I    247	  I 11.74  V 11.34  L 10.93  K  7.69  P  7.69  E  6.88  R  5.67  F  5.26  A  4.86  Q  4.45  T  4.05  D  4.05  S  2.83  Y  2.83  W  2.43  N  2.43  H  1.62  G  1.21  M  0.81  C  0.81  X  0.40
    with open("output_freq.dat", "w") as f:
        for key in amino_acid_counts:
            f.write(f"DATA: {key}")
            #total = sum(total_counts[key].values())
            for amino_acid in sorted(amino_acid_counts[key], key=lambda aa: amino_acid_counts[key][aa], reverse=True):
                count = amino_acid_counts[key][amino_acid]
                total = total_counts[key][amino_acid]
                freq = count / total * 100
                f.write(f",{amino_acid} {freq:.2f} {count}")    
            f.write("\n")
            
    with open("output_freq_select.dat", "w") as f:
        for key in amino_acid_counts:
            
            key_amino_acid = key.split("_")[1]
            if key_amino_acid in amino_acid_counts[key] and key_amino_acid in total_counts[key]:
                key_freq = amino_acid_counts[key][key_amino_acid] / total_counts[key][key_amino_acid] * 100
            else:
                key_freq = 0
                
          

            #total = sum(total_counts[key].values())
            first = True
            for amino_acid in sorted(amino_acid_counts[key], key=lambda aa: amino_acid_counts[key][aa], reverse=True):
                count = amino_acid_counts[key][amino_acid]
                total = total_counts[key][amino_acid]
                freq = count / total * 100
                
                if freq - key_freq >= 10 and count >= 50:
                    
                    if first:
                        first = False
                        f.write(f"SELECT: {key}")
                        f.write(f" {key_freq:.2f}")  
                        
                    f.write(f",{amino_acid} {freq:.2f} {count}")
            if not first:
                f.write("\n")
            
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python aaScanner.py <pdb_file> <databases_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    # Get all residue numbers for chain A into a list
    residue_numbers = []
    with open(pdb_file, "r") as pdb:
        for line in pdb:
            if line.startswith("ATOM") and line[21] == "A":  # Chain identifier is at position 22 (index 21)
                resi = int(line[22:26].strip())  # Residue number is at position 23-26 (index 22-26)
                if resi not in residue_numbers:
                    residue_numbers.append(resi)
    
    # Sort the residue numbers
    residue_numbers.sort()
    
    # Make selection statements like "chain A and resi X-Y"
    selection_statements = []
    start = residue_numbers[0]
    end = start
    
    step = max(1, len(residue_numbers) // 10)
    for i in range(0, len(residue_numbers), step):
        end = residue_numbers[min(i + step - 1, len(residue_numbers) - 1)]
        selection_statements.append(f"chain A and resi {residue_numbers[i]}-{end}")
    
    
    # Print selection statements
    for statement in selection_statements:
        print(statement)
    
    
    databases_file = sys.argv[2]

    with open(databases_file, "r") as file:
        databases = [line.strip() for line in file if line.strip()]  # Read database paths from the file

    main(pdb_file, databases, selection_statements)
