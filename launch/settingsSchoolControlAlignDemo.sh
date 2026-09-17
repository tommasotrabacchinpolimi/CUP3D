#!/bin/bash
# Hybrid Wang school demo: force surge + act({b}) yaw.
# Mesh: extent=1, bpd=2/2/1, Lmax=6 → minH≈0.00195 (~100 cells/L; safe for amp=1 act).

NNODE=${NNODE:-1}
LAMBDA=${LAMBDA:-1e9}
BPDX=${BPDX:-2}
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-1}
LEVELS=${LEVELS:-6}
LEVEL0=${LEVEL0:-3}
CFL=${CFL:-0.3}
PT=${PT:-1e-5}
PTR=${PTR:-1e-3}
NU=${NU:-0.00004}
TEND=${TEND:-3.0}
ZPOS=${ZPOS:-0.25}
EXTENT=${EXTENT:-1.0}

FACTORY=
# Approach+align from rest:
#   xvel=0, sep≈0.28 (~1.45 dAtt), Δφ=12°, γ_R=0.05, amp=1.0
FACTORY+="StefanFish L=0.2 T=0.5 amplitudeFactor=1.0 xvel=0 xpos=0.40 ypos=0.49 zpos=${ZPOS} planarAngle=0 bFixToPlanar=1 bFixFrameOfRef_x=1 CorrectPosition=0 CorrectPositionZ=0 CorrectRoll=0 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=0.5 amplitudeFactor=1.0 xvel=0 xpos=0.67 ypos=0.53 zpos=${ZPOS} planarAngle=12 bFixToPlanar=1 CorrectPosition=0 CorrectPositionZ=0 CorrectRoll=0 heightProfile=danio widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -poissonSolver iterative"
OPTIONS+=" -extent ${EXTENT} -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -tdump 0 -tend ${TEND}"
OPTIONS+=" -CFL ${CFL} -lambda ${LAMBDA} -nu ${NU}"
OPTIONS+=" -levelMax ${LEVELS} -levelStart ${LEVEL0} -levelMaxVorticity ${LEVELS} -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -schoolControl 1"
OPTIONS+=" -schoolK 1 -schoolGammaAtt 0.04 -schoolGammaAli 0.20 -schoolGammaR 0.05"
OPTIONS+=" -schoolDecide 0.1 -schoolV0Cm 14"
OPTIONS+=" -schoolKpPhi 8.0 -schoolKdPhi 1.5 -schoolKpV 6.0"
OPTIONS+=" -schoolBMax 2.0 -schoolAMax 4.0"
OPTIONS+=" -schoolSeed 71 -schoolLogEvery 10"
OPTIONS+=" -verbose 1"
