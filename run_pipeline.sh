#!/bin/bash

#SBATCH -J nf-scrna-seq
#SBATCH -p general
#SBATCH -o run_pipeline_%j.txt
#SBATCH -e run_pipeline_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

echo "Starting Nextflow pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

base=/N/project/Krolab/Siddharth/Pipelines/scrna-seq
cd $base

# Create output directory variable (fixed case)
OUTPUT_DIR=${base}/results                   # Fixed: Changed from lowercase 'output'
mkdir -p ${OUTPUT_DIR}/reports               # Create directory for reports

echo "Base directory: $base"
echo "Output directory: $OUTPUT_DIR"
export NXF_WORK="${base}/work"

# Load modules
module load fastqc
module load sra-toolkit
module load cellranger
module load java/17.0.7
module load conda

echo "Modules loaded successfully"

# Activate conda environment
echo "Activating conda environment..."
conda activate conda_envs/env/renv

echo "Starting pipeline execution..."

nextflow run main.nf \
    --samplesheet samples.csv \
    -resume \
    -profile slurm \
    -with-timeline ${OUTPUT_DIR}/pipeline_info/execution_timeline_${SLURM_JOB_ID}.html \
    -with-report ${OUTPUT_DIR}/pipeline_info/execution_report_${SLURM_JOB_ID}.html \
    -with-trace ${OUTPUT_DIR}/pipeline_info/execution_trace_${SLURM_JOB_ID}.txt \
    -with-dag ${OUTPUT_DIR}/pipeline_info/pipeline_dag_${SLURM_JOB_ID}.html

# Capture exit status
EXIT_STATUS=$?

echo ""
echo "Pipeline finished at $(date)"
echo "Exit status: $EXIT_STATUS"

if [ $EXIT_STATUS -eq 0 ]; then
    echo "SUCCESS: Pipeline completed successfully!"
    echo "Results available at: $OUTPUT_DIR"
    echo "Reports available at: ${OUTPUT_DIR}/reports/"
else
    echo "FAILED: Pipeline failed with exit status $EXIT_STATUS"
    echo "Check logs for details"
fi

exit $EXIT_STATUS