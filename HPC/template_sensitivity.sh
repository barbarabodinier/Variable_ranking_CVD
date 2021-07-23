#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N sens
#PBS -q med-bio
#PBS -J 1-100

cd /work/bbodinie/UK_Biobank/CVD_biomarkers/Scripts
module load anaconda3/personal

Rscript stability_sensitivity.R {data_input} {outcome_input} {model_id_input} {gender_input} $PBS_ARRAY_INDEX
