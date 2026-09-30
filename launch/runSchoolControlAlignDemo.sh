#!/bin/bash
# Scheduler-agnostic runner for the hybrid school align demo.
# Submit via:
#   sbatch school_align11.slurm          # Slurm
#   qsub  school_align11.pbs             # PBS (Aurora)
# Or inside an interactive allocation. Do NOT run on a login node unless FORCE_LOGIN=1.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SETTINGS="${ROOT}/launch/settingsSchoolControlAlignDemo.sh"
BIN="${ROOT}/build/cubismup3d_simulation"
RUNDIR="${ROOT}/runs/school_control_align11"

if [[ -z "${FORCE_LOGIN:-}" ]] && [[ -z "${PBS_JOBID:-}${SLURM_JOB_ID:-}${COBALT_JOBID:-}" ]]; then
  echo "Refusing to run CFD on what looks like a login node."
  echo "Submit school_align11.slurm (or .pbs), or: FORCE_LOGIN=1 $0"
  exit 1
fi

if [[ ! -x "${BIN}" ]]; then
  echo "Missing binary: ${BIN}"
  echo "Build first (cmake + make) on this cluster."
  exit 1
fi

source "${SETTINGS}"
mkdir -p "${RUNDIR}"
rm -f "${RUNDIR}"/field.restart \
      "${RUNDIR}"/*.restart \
      "${RUNDIR}"/*.h5 \
      "${RUNDIR}"/*.xmf \
      "${RUNDIR}"/school_control.log 2>/dev/null || true

cp "${SETTINGS}" "${RUNDIR}/settings.sh"
cp "${BIN}" "${RUNDIR}/simulation"
cd "${RUNDIR}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
NP="${NP:-${SLURM_NTASKS:-1}}"
TEND="${TEND:-2.5}"

echo "NP=${NP}  OMP_NUM_THREADS=${OMP_NUM_THREADS}  tend=${TEND}  (align11 hybrid F+act@0.5T, Lmax=6)"
echo "${OPTIONS}" > settings.txt

# Prefer srun under Slurm; otherwise mpirun (PBS / interactive).
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v srun >/dev/null 2>&1; then
  echo "launcher: srun -n ${NP}"
  srun -n "${NP}" --cpus-per-task="${OMP_NUM_THREADS}" \
    ./simulation ${OPTIONS} -factory-content "${FACTORY}"
else
  echo "launcher: mpirun -np ${NP}"
  mpirun -np "${NP}" ./simulation ${OPTIONS} -factory-content "${FACTORY}"
fi

echo "Log: ${RUNDIR}/school_control.log"
tail -n 30 school_control.log || true
