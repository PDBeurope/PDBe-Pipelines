#!/usr/bin/env python3
"""
Dummy script to copy file contents from input to output.
Used as a placeholder for Nextflow pipeline steps.

Usage:
    python copy_file.py --input input_file.txt --output output_file.txt
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Copy file contents from input to output"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input file path"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file path"
    )

    args = parser.parse_args()

    try:
        with open(args.input, 'r') as input_file:
            contents = input_file.read()
        
        with open(args.output, 'w') as output_file:
            output_file.write(contents)
        
        print(f"Successfully copied contents from {args.input} to {args.output}")
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error copying file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
