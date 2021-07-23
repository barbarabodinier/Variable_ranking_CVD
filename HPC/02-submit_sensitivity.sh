array=(1 3 7)

for j in $(seq 1 1 3)
do
echo $j
echo ${array[$j]}

sed "s/{data_input}/updated/g" template_sensitivity.sh > run0.sh
sed "s/{model_id_input}/${array[$j]}/g" run0.sh > run1.sh
sed "s/{outcome_input}/CVD/g" run1.sh > run2.sh
sed "s/{gender_input}/0/g" run2.sh > run.sh
qsub run.sh

sed "s/{data_input}/updated/g" template_sensitivity.sh > run0.sh
sed "s/{model_id_input}/${array[$j]}/g" run0.sh > run1.sh
sed "s/{outcome_input}/CVD/g" run1.sh > run2.sh
sed "s/{gender_input}/1/g" run2.sh > run.sh
qsub run.sh

done
