#!/bin/bash
# Free planar two-fish + Wang–PD school control.
# Social law / kick noise / cruise / decide period follow Wang et al. 2022.
# Only the outer PD (Kp/Kd → act) is extra (paper has no Stefan act layer).

NNODE=${NNODE:-1}
LAMBDA=${LAMBDA:-1e9}
# Bigger xy mesh, thin z: bpdz=1 → extentz = extent * bpdz/max(bpd).
# bpdx=bpdy=4, bpdz=1, extent=2 → domain 2.0 × 2.0 × 0.50
# (8/8/1 + levelMax=6 blew up to 262k blocks; this is ~2× finer h than the
# old 2/2/1 L=5 run while keeping z half as thick and AMR start moderate.)
BPDX=${BPDX:-4}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-1}
LEVELS=${LEVELS:-6}
LEVEL0=${LEVEL0:-3}
CFL=${CFL:-0.3}
PT=${PT:-1e-5}
PTR=${PTR:-1e-3}
NU=${NU:-0.00004}
TEND=${TEND:-5.0}
ZPOS=${ZPOS:-0.25}

FACTORY=
# bFixToPlanar: z + roll/pitch locked, yaw free.
# T=0.5 + amplitudeFactor=1.5: stronger Stefan plant. xvel=0: no IC speed boost;
# free dynamics from rest (Obstacle.cpp only applies xvel as initial transVel).
FACTORY+="StefanFish L=0.2 T=0.5 amplitudeFactor=1.5 xvel=0 xpos=0.50 ypos=0.50 zpos=${ZPOS} planarAngle=0 bFixToPlanar=1 bFixFrameOfRef_x=1 CorrectPosition=0 CorrectPositionZ=0 CorrectRoll=0 heightProfile=danio widthProfile=stefan
"
FACTORY+="StefanFish L=0.2 T=0.5 amplitudeFactor=1.5 xvel=0 xpos=0.80 ypos=0.55 zpos=${ZPOS} planarAngle=15 bFixToPlanar=1 CorrectPosition=0 CorrectPositionZ=0 CorrectRoll=0 heightProfile=danio widthProfile=stefan
"

OPTIONS=
OPTIONS+=" -poissonSolver iterative"
OPTIONS+=" -extent 2.0 -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -tdump 0 -tend ${TEND}"
OPTIONS+=" -CFL ${CFL} -lambda ${LAMBDA} -nu ${NU}"
OPTIONS+=" -levelMax ${LEVELS} -levelStart ${LEVEL0} -levelMaxVorticity ${LEVELS} -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -poissonTol ${PT} -poissonTolRel ${PTR} "
OPTIONS+=" -schoolControl 1"
# Wang Fig. 5A schooling + paper kick noise / mean inter-kick ~0.5s / v0=14 cm/s
OPTIONS+=" -schoolK 1 -schoolGammaAtt 0.04 -schoolGammaAli 0.20 -schoolGammaR 0.20"
OPTIONS+=" -schoolDecide 0.5 -schoolV0Cm 14"
# PD only (not in paper) — maps setpoints → act({b,a})
# KpV: 0.08 → 0.40 (5×) → 3.2 (8× again) so a pegs near -aMax when U ≪ v_des
OPTIONS+=" -schoolKpPhi 1.2 -schoolKdPhi 0.25 -schoolKpV 3.2"
OPTIONS+=" -schoolBMax 1.5 -schoolAMax 0.5"
OPTIONS+=" -schoolSeed 71 -schoolLogEvery 20"
OPTIONS+=" -verbose 1"
