#!/bin/bash
for i in $@; do
    mom=$(printf "%.2f" ${i});
    cp -f perdaix_pions.mac perdaix_pions_${mom}_GeV.mac;
    sed -i "s/MOMENTUM/${mom}/g" perdaix_pions_${mom}_GeV.mac;
done
    
    