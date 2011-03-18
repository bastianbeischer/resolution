#!/bin/bash
for i in $@; do
    mom=$(printf "%.2f" ${i});
    cp -f perdaix_protons.mac perdaix_protons_${mom}_GeV.mac;
    sed -i "s/MOMENTUM/${mom}/g" perdaix_protons_${mom}_GeV.mac;
done
    
    