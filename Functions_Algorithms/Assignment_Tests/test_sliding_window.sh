#!/usr/bin/env bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo "Running Python unit tests first:"
python3 test_sliding_window.py

echo "Custom bash-based tests:"
python3 sliding_window.py &>/dev/null
if [ $? -eq 1 ]; then
    echo -e "${GREEN}PASSED:${NC} No arguments should cause an error"
    else
    echo -e "${RED}FAILED:${NC} No arguments should cause an error"
fi

output_example=$(python3 sliding_window.py 13 ATCTCTAGCTGATGCTGATGCTGATGCTGATGAATATCGCGCGATATATATCGATCGATCGATGGCGATCGATCGAGCTAGCGATCTCTCCTGATGCGTA 2>/dev/null)
output_correct=$(cat sliding_window_testout)
if [ "$output_example" == "$output_correct" ]; then
    echo -e "${GREEN}PASSED:${NC} Random nucleotide string with 13-mers"
    else
    echo -e "${RED}FAILED:${NC} Random nucleotide string with 13-mers"
    printf "\tSee expected output in sliding_window_testout, in this directory\n"
fi
