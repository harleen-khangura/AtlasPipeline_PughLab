import sys
import re

def convert_bed(input_bed, output_bed):
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            if line.strip():  # Skip empty lines
                fields = line.split('\t')
                
                # Extract relevant fields
                chromosome_info = fields[0]
                start_offset = int(fields[1])  # Start position offset as integer
                end_offset = int(fields[2])  # End position offset as integer
                strand = fields[-1]  # Assuming strand is the last field

                # Extract the base position from the chromosome information
                match = re.match(r'([^:]+):(\d+)-(\d+)', chromosome_info)
                if match:
                    chromosome = match.group(1)
                    base_start = int(match.group(2))  # Base start position from chromosome info
                    base_end = int(match.group(3))  # Base end position from chromosome info

                    # Calculate new start and end positions based on strand
                    if strand == '+':
                        new_start = base_start + start_offset
                        new_end = base_start + end_offset
                    elif strand == '-':
                        new_start = base_end - end_offset
                        new_end = base_end - start_offset
                    else:
                        print(f"Warning: Unexpected strand value: {strand}")
                        continue
                else:
                    print(f"Warning: Could not parse chromosome info: {chromosome_info}")
                    continue
                
                # Create new line for output
                new_line = f"{chromosome}\t{new_start}\t{new_end}\t.\t.\t{strand}\n"
                outfile.write(new_line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_bed.py <input_bed_file> <output_bed_file>")
        sys.exit(1)

    input_bed_file = sys.argv[1]
    output_bed_file = sys.argv[2]

    convert_bed(input_bed_file, output_bed_file)
    print(f"Converted BED file saved as: {output_bed_file}")
