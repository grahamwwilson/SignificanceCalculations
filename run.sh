#!/bin/sh

# Choose order of log file
#                       PULL              NOBS    NEXP
./testPoissonArgs 170  5.213856  5.227546  13      3.38383
./testPoissonArgs 330  4.689327  4.733014  35     16.04270
./testPoissonArgs 30  3.857803  3.883885  38      20.44050
./testPoissonArgs 70 -3.586047 -3.912632  787    904.684
./testPoissonArgs 18  3.496100  3.708803  113     79.85710

./testPoissonArgs 60 -3.467514 -3.545179  118    163.304
./testPoissonArgs 40 -3.224744 -3.311679  27      50.5442
./testPoissonArgs 10  3.146387  3.193373  21      10.6021
./testPoissonArgs 10  3.080125  3.270015  159    122.768
./testPoissonArgs 10  2.990918  3.011731  10       3.98662

exit
