### Split a file with two columns using the values in the second column as the unique id for the split
### use this to create the input pop files to check missing site data per population

# Read data from the input text file with two columns
with open('wgs_all.txt', 'r') as file:
    lines = file.readlines()

# Create a dictionary to store data based on unique values in the second column
data_dict = {}
for line in lines:
    # Split the line into columns
    columns = line.strip().split('\t')  # Assuming columns are separated by tabs, adjust if necessary
    
    # Extract values from the columns
    first_column = columns[0]
    second_column = columns[1].strip()
    
    # Check if the second column value is already a key in the dictionary
    if second_column in data_dict:
        # If it is, append the first column value to the existing list of values for that key
        data_dict[second_column].append((first_column, second_column))
    else:
        # If it is not, create a new key-value pair with the second column value as the key
        # and a list containing the first column value as the initial value
        data_dict[second_column] = [(first_column, second_column)]

# Write data to separate text files named after the unique values in the second column
for key, values in data_dict.items():
    output_file_name = f'{key}.txt'
    with open(output_file_name, 'w') as output_file:
        for value_pair in values:
            output_file.write(f"{value_pair[0]}\t{value_pair[1]}\n")
