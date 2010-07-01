#!/bin/bash

for i in $(find . -name "*000*" | sort); do
    stem=$(echo $i | sed 's/\(.*\)\(_000\)\(\.root\)/\1/');
    hadd ${stem}.root ${stem}_00[0-9].root;
done;