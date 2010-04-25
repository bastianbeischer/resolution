#!/bin/bash

NTHREADS=8;

E=( 010 020 050 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 );

for (( i=0; i<${#E[*]}; )); do
    N=`ps axu | grep resolution | grep -v grep | wc -l`;
    if [[ N -lt NTHREADS ]]; then
        filestem=pebs01_${E[i]}_GeV;
        resolution mac/${filestem}.mac >& out/${filestem}.out&
        i=$(expr $i + 1);
    fi;
    sleep 1;
done
