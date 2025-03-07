#!/usr/bin/env python

import sys

def calculate_mean(input_file, output_file):
    try:
        # Read numbers from input file
        numbers = []
        with open(input_file, 'r') as f:
            for line in f:
                # Strip whitespace and convert to float
                num = float(line.strip())
                numbers.append(num)
        
        # Check if list is empty
        if not numbers:
            raise ValueError("Input file is empty")
        
        # Calculate mean
        mean = sum(numbers) / len(numbers)
        
        # Round to 2 decimal places
        mean_rounded = round(mean, 3)
        
        # Write result to output file
        with open(output_file, 'w') as f:
            f.write(str(mean_rounded))
            
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
    except ValueError as e:
        if str(e) == "Input file is empty":
            print("Error: Input file is empty")
        else:
            print("Error: Invalid number format in input file")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    # Check if correct number of arguments provided
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)
    
    # Get input and output filenames from command line arguments
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
    # Calculate mean and save result
    calculate_mean(input_filename, output_filename)