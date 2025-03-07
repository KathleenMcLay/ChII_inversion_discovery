### Submit multiple GNU parallel jobs 

import subprocess

# Open the text file containing the list of samples and pops 
file_path = '/home/564/km6006/Scripts/inversion_paper/wgs_all.txt'
with open(file_path, 'r') as file:
    # Read lines from the file
    lines = file.readlines()

# Create a list of lists comprised of the first element from each row (so the values in column 1)
result = [line.strip().split('\t')[0] for line in lines]
# Create sublists of 10 commands
sublists = [result[i:i + 10] for i in range(0, len(result), 10)]

job_num = 0

for i in sublists:
    job_num = job_num + 1

    job_string = f"""#!/bin/bash -l
#PBS -N data_prep_{job_num}
#PBS -l walltime=48:00:00
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l jobfs=400GB
#PBS -o pbs_job_output_{job_num}.out
#PBS -l storage=scratch/ht96+gdata/ht96

module load parallel
parallel --jobs 5 /home/564/km6006/Scripts/inversion_paper/2_data_prep.sh ::: {' '.join(i)} &>> /home/564/km6006/Scripts/inversion_paper/2_data_prep_submit/2_data_prep_out_{job_num}.txt  

"""

    # Print the job string instead of running it to test
    print(job_string)

    # For testing, you might want to comment out the next lines
    # # Write the job string to a temporary script file
    script_file_path = 'temp_script.sh'
    with open(script_file_path, 'w') as script_file:
        script_file.write(job_string)

    # Run the script using subprocess
    subprocess.run(['qsub', script_file_path])

    # Clean up the temporary script file
    subprocess.run(['rm', script_file_path])
