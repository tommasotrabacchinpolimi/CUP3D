#!/bin/bash
# Free planar two-fish Wang–PD school control.
# Do NOT run on a login node unless FORCE_LOGIN=1.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SETTINGS="${ROOT}/launch/settingsSchoolControlTwoFish.sh"
BIN="${ROOT}/build/cubismup3d_simulation"
RUNDIR="${ROOT}/runs/school_control_free"

if [[ -z "${FORCE_LOGIN:-}" ]] && [[ -z "${PBS_JOBID:-}${SLURM_JOB_ID:-}${COBALT_JOBID:-}" ]]; then
  echo "Refusing to run CFD on what looks like a login node."
  echo "Submit this script inside a job, or: FORCE_LOGIN=1 $0"
  exit 1
fi

source "${SETTINGS}"
mkdir -p "${RUNDIR}"
# Stale field.restart / HDF5 dumps make init deserialize a broken mesh → GSL singular.
rm -f "${RUNDIR}"/field.restart \
      "${RUNDIR}"/*.restart \
      "${RUNDIR}"/*.h5 \
      "${RUNDIR}"/*.xmf \
      "${RUNDIR}"/school_control.log 2>/dev/null || true

cp "${SETTINGS}" "${RUNDIR}/settings.sh"
cp "${BIN}" "${RUNDIR}/simulation"
cd "${RUNDIR}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}  tend=${TEND}  (fresh rundir, no restart)"
echo "${OPTIONS}" > settings.txt
mpirun -np 1 ./simulation ${OPTIONS} -factory-content "${FACTORY}"
echo "Log: ${RUNDIR}/school_control.log"
tail -n 20 school_control.log || true
