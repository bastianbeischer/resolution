#!/usr/bin/perl -w

use strict;

#
# EDIT firstRunNumber
# EDIT enter VALUEs
# EDIT macro template to place values
# RUN script
#

my $firstRunNumber = 1000;
my $numberOfEventsPerRun = 1000;
my $maximumEnergy = 100; # GeV

my $condor_dir = "/home/home4/institut_1b/beischer/src/geant4/resolution/condor";
my $executable = "/home/home4/institut_1b/beischer/src/geant4/resolution/bin/moontrace";

################################################################################################################

my $currentRun = $firstRunNumber;

for (my $energy = 1; $energy <= $maximumEnergy; ++$energy) {
    print "Run $currentRun:\n";
    my $condorfile = &make_condor_file($currentRun);
    print "$condorfile\n";
    my $macrofile = &make_macro_file($currentRun, $energy, $numberOfEventsPerRun);
    print "$macrofile\n";

    unlink("/home/home4/institut_1b/beischer/src/geant4/resolution/results/res_${currentRun}.root");
    system "condor_submit", "$condorfile";
    ++$currentRun;
}

sub make_macro_file{
  # make a macro file, first argument is run number, second is current value (to be placed in macro template below)

  my ($currentRun, $energy, $nEvents) = @_;

  open MACROFILE, ">${condor_dir}/mac/res_${currentRunNumber}.mac" or die "Error: Cannot make macro file: $!";

  print MACROFILE <<EOF;

/control/verbose 1
/run/verbose 0
/tracking/verbose 0
/RES/Fit/Verbose 0

/RES/Det/AddModule 0. 0. 45. cm
/RES/Det/AddModule 0. 0. 40. cm
/RES/Det/AddModule 0. 0. 35. cm
/RES/Det/AddModule 0. 0. 30. cm
/RES/Det/AddModule 0. 0. 25. cm
/RES/Det/AddModule 0. 0. 20. cm
/RES/Det/AddModule 0. 0. 15. cm
/RES/Det/AddModule 0. 0. 10. cm
/RES/Det/AddModule 0. 0. 5. cm
/RES/Det/AddModule 0. 0. 0. cm
/RES/Det/AddModule 0. 0. -5. cm
/RES/Det/AddModule 0. 0. -10. cm
/RES/Det/AddModule 0. 0. -15. cm
/RES/Det/AddModule 0. 0. -20. cm
/RES/Det/AddModule 0. 0. -25. cm
/RES/Det/AddModule 0. 0. -30. cm
/RES/Det/AddModule 0. 0. -35. cm
/RES/Det/AddModule 0. 0. -40. cm
/RES/Det/AddModule 0. 0. -45. cm
/RES/Det/ModuleGap 2.4 cm
/RES/Det/ModuleWidth 99.5 cm
/RES/Det/ModuleLength 99.5 cm

/RES/Gun/Energy $energy

#/RES/Field/SetInhomFieldFrom tables/perdaix_07_jul_2009.table
#/RES/Field/SetDummyField 0.3 0.0 0.0 tesla
/RES/Field/SetUniformField 0.3 0.0 0.0 tesla

/RES/Data/OverWriteFile true
/RES/Data/SetFileName /home/home4/institut_1b/beischer/src/geant4/resolution/results/res_$currentRun.root
/RES/Run/StoreResults

/run/initialize

# /vis/scene/create
# /vis/scene/add/volume world
# /vis/scene/add/trajectories

# /vis/open OGLIX

# /vis/scene/endOfEventAction accumulate
# /vis/viewer/set/viewpointThetaPhi 90 0
# /vis/viewer/set/upVector 0 0 1
# /vis/viewer/zoom 1.5

/RES/Run/Generate $nEvents
/RES/Run/Reconstruct

EOF

  close MACROFILE;

  # return value: macro file name
  "${condor_dir}/mac/res_${currentRunNumber}.mac";

}



sub make_condor_file{
  # make a condor file, first ( and only ) argument is run number

  my $ith = shift;
  my $currentRunNumber = $ith;

  open CONDORFILE, ">${condor_dir}/condor/res_${currentRunNumber}.condor" or die "Error: Cannot make condor file: $!";

  print CONDORFILE <<EOF;
universe        = vanilla
executable      = $executable
arguments	= ${condor_dir}/mac/res_$currentRunNumber.mac
initialdir      = ${condor_dir}
EOF

  print CONDORFILE <<'EOF';
getenv          = true
requirements    = ( OpSys == "LINUX" )
rank            = KFlops
EOF

  print CONDORFILE <<EOF;
error           = ${condor_dir}/ERR/res_${currentRunNumber}.ERR
output          = ${condor_dir}/STD/res_${currentRunNumber}.STD
log             = ${condor_dir}/LOG/res_${currentRunNumber}.LOG
EOF

  print CONDORFILE <<'EOF';
notification    = error
notify_user	= beischer@physik.rwth-aachen.de
coresize        = 1000
Queue
EOF

  close CONDORFILE;

  # return value: condor file name
  "${condor_dir}/condor/res_${currentRunNumber}.condor";

}
