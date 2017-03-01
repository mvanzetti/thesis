#!/usr/bin/env bash

echo Automated overlap and reldist analysis

cd ..

echo ------ Computing Relative Distances
#./analyze.py ENCODE FANTOM --type "tissue" --analysis reldist
#echo
#./analyze.py ENCODE FANTOM --type "stem cell" --analysis reldist
#echo
./analyze.py ENCODE FANTOM --type "primary cell" --analysis reldist
#echo
#./analyze.py ENCODE FANTOM --type "immortalized cell line" --analysis reldist
#echo
#./analyze.py ENCODE FANTOM --type "in vitro differentiated cells" --analysis reldist
#echo
#./analyze.py ENCODE FANTOM --type "induced pluripotent stem cell line" --analysis reldist

echo
echo
echo ------- Computing Overlap Analysis
#./analyze.py ENCODE FANTOM --intervals 20 --samples 50 --type "tissue" --analysis overlap
#echo
#./analyze.py ENCODE FANTOM --intervals 20 --samples 50 --type "stem cell" --analysis overlap
echo
./analyze.py ENCODE FANTOM --intervals 20 --samples 50 --type "primary cell" --analysis overlap
echo
./analyze.py ENCODE FANTOM --intervals 20 --samples 50 --type "immortalized cell line" --analysis overlap
echo
./analyze.py ENCODE FANTOM --intervals 20 --samples 50 --type "in vitro differentiated cells" --analysis overlap
echo
./analyze.py ENCODE FANTOM --intervals 20 --samples 50 --type "induced pluripotent stem cell line" --analysis overlap