PROGRAM="rsync"
ARGUMENTS="-avz --exclude=CVS"
REMOTEPATH="beischer@portal.physik.rwth-aachen.de:~/src/geant4/resolution"

${PROGRAM} ${ARGUMENTS} ${REMOTEPATH}/results/*.root results/
${PROGRAM} ${ARGUMENTS} ${REMOTEPATH}/condor/condor.pl condor/
${PROGRAM} ${ARGUMENTS} ${REMOTEPATH}/condor/condor condor/
${PROGRAM} ${ARGUMENTS} ${REMOTEPATH}/condor/ERR condor/
${PROGRAM} ${ARGUMENTS} ${REMOTEPATH}/condor/LOG condor/
${PROGRAM} ${ARGUMENTS} ${REMOTEPATH}/condor/mac condor/
${PROGRAM} ${ARGUMENTS} ${REMOTEPATH}/condor/STD condor/