#!/bin/bash
for i in $@; do
    mom=$(printf "%.2f" ${i});
    cp -f perdaix_helium.mac perdaix_helium_${mom}_GeV.mac;
    sed -i "s/MOMENTUM/${mom}/g" perdaix_helium_${mom}_GeV.mac;
done
    
    