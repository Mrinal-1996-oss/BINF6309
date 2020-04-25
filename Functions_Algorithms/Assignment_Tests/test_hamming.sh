#!/usr/bin/env bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo "Running Python unit tests first:"
python3 test_hamming.py

echo "Custom bash-based tests:"

python3 hamming.py &>/dev/null
if [ $? -eq 1 ]; then
    echo -e "${GREEN}PASSED:${NC} No arguments should cause an error"
    else
    echo -e "${RED}FAILED:${NC} No arguments should cause an error"
fi

output_example=$(python3 hamming.py ATCGATCGATCGGGCTA AAAAAAAAAAAAAAAAA 2>/dev/null)
output_example_compare=$(echo $output_example | sed -r 's/\s+//g')
output_correct="ATCGATCGATCGGGCTA\tAAAAAAAAAAAAAAAAA\t13" 
output_correct_compare=$(printf $output_correct | sed -r 's/\s+//g')
if [ $output_example_compare == $output_correct_compare ]; then
    echo -e "${GREEN}PASSED:${NC} Random 17-mer comparison"
    else
    echo -e "${RED}FAILED:${NC} Random 17-mer compraison"
    printf "\tActual: ${output_example} Expected: ${output_correct}\n"
    echo $output_example_compare
    echo $output_correct_compare
fi
