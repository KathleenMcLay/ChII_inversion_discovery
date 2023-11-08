### Split a text file into sub files based on number of rows
### Use this to split the GNU parallel txt files and run multiple parallel jobs

# Define the input and output file names
input_file = 'vc_file_list.txt'
output_file_prefix = 'variant_calling_list_'

# Function to split text into chunks and write to files
def split_and_write(input_file, chunk_size):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Calculate the number of chunks
    num_chunks = len(lines) // chunk_size + (len(lines) % chunk_size != 0)

    # Split the lines into chunks and write to separate files
    for i in range(num_chunks):
        chunk_start = i * chunk_size
        chunk_end = (i + 1) * chunk_size
        chunk = lines[chunk_start:chunk_end]

        # Write the chunk to a new file
        output_file_name = f'{output_file_prefix}{i + 1}.txt'
        with open(output_file_name, 'w') as output_file:
            output_file.writelines(chunk)

# Call the function with the desired chunk size
split_and_write(input_file, 10)