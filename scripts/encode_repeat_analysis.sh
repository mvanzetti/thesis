#!/usr/bin/env bash

echo ENCODE and RepeatMasker - Automated overlap and reldist analysis

cd ..

echo ------ Computing Relative Distances
./analyze.py ENCODE RepeatMasker --type "tissue" --repeat_class "SINE/MIR" --analysis reldist
echo
./analyze.py ENCODE RepeatMasker --type "stem cell" --repeat_class "SINE/MIR" --analysis reldist
echo
#./analyze.py ENCODE RepeatMasker --type "primary cell" --repeat_class "SINE/Alu" --analysis reldist
#echo
#./analyze.py ENCODE RepeatMasker --type "immortalized cell line"  --repeat_class "SINE/Alu" --analysis reldist
#echo
#./analyze.py ENCODE RepeatMasker --type "in vitro differentiated cells" --repeat_class "SINE/Alu" --analysis reldist
#echo
#./analyze.py ENCODE RepeatMasker --type "induced pluripotent stem cell line" --repeat_class "SINE/Alu" --analysis reldist



echo
echo
echo ------- Computing Overlap Analysis
./analyze.py ENCODE RepeatMasker --intervals 20 --samples 50 --type "tissue" --repeat_class "SINE/MIR" --analysis overlap
echo
./analyze.py ENCODE RepeatMasker --intervals 20 --samples 50 --type "stem cell" --repeat_class "SINE/MIR" --analysis overlap
echo
#./analyze.py ENCODE RepeatMasker --intervals 20 --samples 50 --type "primary cell" --repeat_class "SINE/Alu" --analysis overlap
#echo
#./analyze.py ENCODE RepeatMasker --intervals 20 --samples 50 --type "immortalized cell line" --repeat_class "SINE/Alu" --analysis overlap
#echo
#./analyze.py ENCODE RepeatMasker --intervals 20 --samples 50 --type "in vitro differentiated cells" --repeat_class "SINE/Alu" --analysis overlap
#echo
#./analyze.py ENCODE RepeatMasker --intervals 20 --samples 50 --type "induced pluripotent stem cell line" --repeat_class "SINE/Alu" --analysis overlap