#!/bin/bash
NNODE=${NNODE:-1}
PSOLVER="iterative" # CPU Poisson solver
#PSOLVER="cuda_iterative" # GPU Poisson solver

# BurstandCoastFish (y(x,t) burst-and-coast midline).
# Kinematics (factory): T=T_bout, lambdaWave and lambdaBout are fractions of L.
# Tail speed ~ 2 pi a c / lambdaWave,  c = lambdaBout * L / T.
#   T=0.3, lambdaWave=0.5, lambdaBout=1.8 => c=1.2, burst duration=0.083
#   full rest while lambdaBout > 1 + lambdaWave (1.8 > 1.5)
# Target U ~ 2 L/s = 0.4, Re = U L / nu = 4000 => nu = 2e-5
# Do not pass -lambda on OPTIONS for the wave (that is penalization).
#
# From launch/:
#   ./launchMac.sh settingsBurstAndCoastFish.sh burstCoast_test
#
# Rebuild first so BurstAndCoastFish.cpp is in bin/simulation.

FACTORY=
FACTORY+="BurstandCoastFish L=0.2 T=0.1 xpos=0.5 ypos=0.5 zpos=0.5 bFixToPlanar=1 bFixFrameOfRef=1 amplitudeFactor=1.0 lambdaWave=0.5 lambdaBout=1.8 bHoldHeading=1 heightProfile=baseline widthProfile=baseline
"
OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -extent 2.0 -bpdx 4 -bpdy 2 -bpdz 2"
OPTIONS+=" -tdump 0.1 -tend 4.0"
OPTIONS+=" -CFL 0.3 -nu 0.00002 -lambda 1e10 "
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4 "
OPTIONS+=" -levelMax 5 -levelStart 3 -levelMaxVorticity 4 -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -restart 0 -checkpointsteps 1000 "
