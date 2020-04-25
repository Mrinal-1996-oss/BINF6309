#!/usr/bin/env bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo "Running Python unit tests first:"
python3 test_fibonacci.py

echo "Custom bash-based tests:"

python3 fibonacci.py &>/dev/null
if [ $? -eq 1 ]; then
    echo -e "${GREEN}PASSED:${NC} No arguments should cause an error"
    else
    echo -e "${RED}FAILED:${NC} No arguments should cause an error"
fi

output_example=$(python3 fibonacci.py 100 100 2>/dev/null)
output_correct="73956552763327946079395804120682457999029562797499247920675380589215558794227968018153512677426569801"
if [ $output_example == $output_correct ]; then
    echo -e "${GREEN}PASSED:${NC} fibonacci 100 100"
    else
    echo -e "${RED}FAILED:${NC} fibonacci 100 100"
    printf "\tActual: ${output_example} Expected: ${output_correct}\n"
    echo $output_example_compare
    echo $output_correct_compare
fi
