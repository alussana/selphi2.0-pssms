#!/usr/bin/env python3

import sys

def main():
    patterns_file = sys.argv[1]
    col = int(sys.argv[2]) - 1

    pattern_list = []
    with open(patterns_file) as patterns:
        for line in patterns:
            pattern_list.append(line.strip())
    
    for line in sys.stdin:
        field = line.split('\t')[col].strip()
        for p in pattern_list:
            if p in field:
                print(line.strip())
                break

if __name__ == '__main__':
    main()
