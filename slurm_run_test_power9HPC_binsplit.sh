#! /bin/bash
#SBATCH -J Fzrefpower
#SBATCH -p main
# test trend with 22-242:22
#SBATCH --array=1-242
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -t 96:00:00
#SBATCH -o logs/Frefpower_%A_%a.out
#SBATCH -e logs/Frefpower_%A_%a.err

# Optional: adjust to your cluster
#module load R || true

# B can be passed via: sbatch --export=B=3000 slurm_ref_power.sh
: "${B:=5000}"
: "${INDIR:=out5}"
: "${TESTOF:=binsplit}"
: "${OUTDIR:=out11}"

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}  B=${B} INDIR=${OUTDIR} TEST=${TESTOF} OUTDIR=${INDIR}"
./run_test_power9HPC.R "${SLURM_ARRAY_TASK_ID}" "${B}" "${OUTDIR}" "${TESTOF}" "${INDIR}"
