#!/usr/bin/env bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

python3 sliding_window_fasta.py &>/dev/null
if [ $? -eq 1 ]; then
    echo -e "${GREEN}PASSED:${NC} No arguments should cause an error"
    else
    echo -e "${RED}FAILED:${NC} No arguments should cause an error"
fi

output_example=$(python3 sliding_window_fasta.py 20 dengue.fasta 2>/dev/null)
output_correct=$(cat sliding_window_fasta_testout)
output_correct_also=$(cat sliding_window_fasta_testout2)
if [ "$output_example" == "$output_correct" ] || [ "$output_example" == "$output_correct_also" ]; then
    echo -e "${GREEN}PASSED:${NC} Dengue virus with 20-mers"
    else
    echo -e "${RED}FAILED:${NC} Dengue virus with 20-mers"
    printf "\tSee expected output in sliding_window_fasta_testout, in this directory\n"
fi
