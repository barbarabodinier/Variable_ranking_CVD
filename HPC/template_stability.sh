#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N lasso
#PBS -q med-bio
#PBS -J 0-1

cd /work/bbodinie/UK_Biobank/CVD_biomarkers/Scripts
module load anaconda3/personal

Rscript stability_selection.R {data_input} {outcome_input} {model_id_input} $PBS_ARRAY_INDEX
