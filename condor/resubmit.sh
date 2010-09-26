#! /bin/bash

for i in $(find STD/ -name "*.STD"); do 
    match=$(tail -n3 $i | head -n1 | grep "Successfully reconstructed")
    if [ -z "${match}" ]; then
        filestem=$(echo $i | sed 's/STD\/\(res_[0-9]*\).STD/\1/');
        condor_submit condor/${filestem}.condor
    fi;
done
