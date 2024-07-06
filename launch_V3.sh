#!/bin/sh

VERSION=V3

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.05 -e 0.015 -s 123450 >test_${VERSION}_0.log
mv out.root out0_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.15 -e 0.015 -s 123451 >test_${VERSION}_0.log
mv out.root out1_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.25 -e 0.015 -s 123452 >test_${VERSION}_0.log
mv out.root out2_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.35 -e 0.015 -s 123453 >test_${VERSION}_0.log
mv out.root out3_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.45 -e 0.015 -s 123454 >test_${VERSION}_0.log
mv out.root out4_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.55 -e 0.015 -s 123455 >test_${VERSION}_0.log
mv out.root out5_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.65 -e 0.015 -s 123456 >test_${VERSION}_0.log
mv out.root out6_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.75 -e 0.015 -s 123457 >test_${VERSION}_0.log
mv out.root out7_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.85 -e 0.015 -s 123458 >test_${VERSION}_0.log
mv out.root out8_${VERSION}.root

time  ./PoissGauss2 -n 1000 -l 40 -u 180   -b 100.95 -e 0.015 -s 123459 >test_${VERSION}_0.log
mv out.root out9_${VERSION}.root

hadd combined_${VERSION}.root out0_${VERSION}.root out1_${VERSION}.root out2_${VERSION}.root out3_${VERSION}.root out4_${VERSION}.root out5_${VERSION}.root out6_${VERSION}.root out7_${VERSION}.root out8_${VERSION}.root out9_${VERSION}.root

exit
