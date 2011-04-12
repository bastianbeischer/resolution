#!/bin/bash
for i in $@; do
    mom=$(printf "%.2f" ${i});
    cp -f perdaix_electrons.mac perdaix_electrons_${mom}_GeV.mac;
    sed -i "s/MOMENTUM/${mom}/g" perdaix_electrons_${mom}_GeV.mac;
done
    
    