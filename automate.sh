#!/bin/sh

N=$1
time ./Parser_WithToys_Toys -t ${N} >Parser_WTT-Standard.log-toy${N}
mv Analyzer.root Analyzer_Toy${N}.root

exit
