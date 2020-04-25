#!/usr/bin/env bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

output_noargs=$(python3 hello_you.py)
output_correct="Hello, you!"
if [ "$output_noargs" == "$output_correct" ]; then
    echo -e "${GREEN}PASSED:${NC} No arguments should return Hello, you!"
    else
    echo -e "${RED}FAILED:${NC} No arguments should return Hello, you!"
    printf "\tExpected: ${output_correct}, Actual: ${output_noargs}\n"
fi

output_example=$(python3 hello_you.py Chesley)
output_correct="Hello, Chesley!"
if [ "$output_example" == "$output_correct" ]; then
    echo -e "${GREEN}PASSED:${NC} 'Chesley' arg should return Hello, Chesley!"
    else
    echo -e "${RED}FAILED:${NC} 'Chesley' arg should return Hello, Chesley!"
    printf "\tExpected: ${output_correct}, Actual: ${output_example}\n"
fi