#!/bin/bash
NNODE=${NNODE:-1}
PSOLVER="iterative" # CPU Poisson solver
#PSOLVER="cuda_iterative" # GPU Poisson solver

# Smoke test for BurstandCoastFish (y(x,t) burst-and-coast midline).
# L=0.2, T=1 => T_bout=1, lambda=0.4 L, lambda_bout=2 L (see BurstAndCoastFish.cpp).
# Re = L^2 / (T nu) = 1000 => nu = 4e-5
#
# From launch/:
#   ./launchMac.sh settingsBurstAndCoastFish.sh burstCoast_test
#
# Rebuild first so BurstAndCoastFish.cpp is in bin/simulation.

FACTORY=
FACTORY+="BurstandCoastFish L=0.2 T=0.5 xpos=0.5 ypos=0.5 zpos=0.5 bFixToPlanar=1 bFixFrameOfRef=1 amplitudeFactor=1.0 heightProfile=baseline widthProfile=baseline
"
OPTIONS=
OPTIONS+=" -poissonSolver ${PSOLVER}"
OPTIONS+=" -extent 2.0 -bpdx 4 -bpdy 2 -bpdz 2"
OPTIONS+=" -tdump 0.1 -tend 4.0"
OPTIONS+=" -CFL 0.3 -nu 0.00004 -lambda 1e10 "
OPTIONS+=" -poissonTol 1e-6 -poissonTolRel 1e-4 "
OPTIONS+=" -levelMax 5 -levelStart 3 -levelMaxVorticity 4 -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -restart 0 -checkpointsteps 1000 "
