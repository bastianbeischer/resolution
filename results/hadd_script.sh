#!/bin/bash

for i in $(find . -regex ".*_electrons_consistency_000\.root$" | sort); do
    stem=$(echo $i | sed 's/\(.*\)\(_000\)\(\.root\)/\1/');
    hadd -f ${stem}.root ${stem}_00[0-9].root;
done;
