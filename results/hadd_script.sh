#!/bin/bash

for i in $(find . -regex ".*_protons_000\.root$" | sort); do
    stem=$(echo $i | sed 's/\(.*\)\(_000\)\(\.root\)/\1/');
    files=${stem}_0[01][0-9].root;
    hadd -f ${stem}.root ${files};
    rm -f ${files};
done;
