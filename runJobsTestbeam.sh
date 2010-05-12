#!/bin/bash

NTHREADS=8;

E=( 00.25 00.50 00.75 01.00 01.50 02.00 02.50 03.00 03.50 04.00 04.50 05.00 05.50 06.00 06.50 07.00 07.50 08.00 08.50 09.00 09.50 10.00 );
#PARTICLES=( electrons pions )

for (( j=0; j<${#E[*]}; )); do
    N=`ps axu | grep resolution | grep -v grep | wc -l`;
    if [[ N -lt NTHREADS ]]; then
        filestem=testbeam_${E[j]}_GeV_electrons_noladder;
        resolution mac/${filestem}.mac >& out/${filestem}.out&
        j=$(expr $j + 1);
    fi;
    sleep 1;
done