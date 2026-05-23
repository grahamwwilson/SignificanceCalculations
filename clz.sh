#!/bin/sh
#
# To compile filename.cpp with ROOT do
# ./cl.sh filename
# Then the executable can be executed using ./filename
#

target=Parser_WithToys
echo 'Compiling with ROOT libraries '${target}.cpp

g++ -g -o ${target} ${target}.cpp `root-config --cflags --glibs`

mv ./Parser_WithToys Zscore.x

exit
