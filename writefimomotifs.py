import argparse
import os

# Define the function to separate lines into GFF files based on MEME value
def separate_meme_lines(input_file, output_dir):
    # Extract the base name of the input file (without extension)
    base_name = os.path.splitext(os.path.basename(input_file))[0]

    # Initialize lists for each MEME category
    meme_1_lines = []
    meme_2_lines = []
    meme_3_lines = []

    # Read the input file
    with open(input_file, 'r') as file:
        for line in file:
            # Check for lines containing Alias=MEME-
            if "Alias=MEME-" in line:
                # Split the line by ';' to find the Alias
                parts = line.split(';')
                for part in parts:
                    if part.startswith("Alias=MEME-"):
                        # Extract the MEME value
                        meme_value = part[len("Alias=MEME-"):].strip()
                        if meme_value == '1':
                            meme_1_lines.append(line)
                        elif meme_value == '2':
                            meme_2_lines.append(line)
                        elif meme_value == '3':
                            meme_3_lines.append(line)
                        break  # No need to check further parts

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Write the separated lines to respective GFF files in the output directory
    with open(os.path.join(output_dir, f'meme_1_{base_name}.gff'), 'w') as meme_1_file:
        meme_1_file.writelines(meme_1_lines)

    with open(os.path.join(output_dir, f'meme_2_{base_name}.gff'), 'w') as meme_2_file:
        meme_2_file.writelines(meme_2_lines)

    with open(os.path.join(output_dir, f'meme_3_{base_name}.gff'), 'w') as meme_3_file:
        meme_3_file.writelines(meme_3_lines)

# Main function to handle command-line arguments
def main():
    parser = argparse.ArgumentParser(description='Separate MEME lines from a GFF file.')
    parser.add_argument('input_file', type=str, help='Path to the input GFF file')
    parser.add_argument('output_dir', type=str, help='Directory to save the output GFF files')
    args = parser.parse_args()

    separate_meme_lines(args.input_file, args.output_dir)

if __name__ == '__main__':
    main()
