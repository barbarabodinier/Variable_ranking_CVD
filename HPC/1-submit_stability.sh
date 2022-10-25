for j in $(seq 1 1 7)
do
echo $j

#sed "s/{model_id_input}/${j}/g" template_stability.sh > run1.sh
#sed "s/{outcome_input}/CAD/g" run1.sh > run.sh
#qsub run.sh

#sed "s/{data_input}/2017/g" template_stability.sh > run0.sh
#sed "s/{model_id_input}/${j}/g" run0.sh > run1.sh
#sed "s/{outcome_input}/CVD/g" run1.sh > run.sh
#qsub run.sh

sed "s/{data_input}/updated/g" template_stability.sh > run0.sh
sed "s/{model_id_input}/${j}/g" run0.sh > run1.sh
sed "s/{outcome_input}/CVD/g" run1.sh > run.sh
qsub run.sh

#sed "s/{data_input}/censored/g" template_stability.sh > run0.sh
#sed "s/{model_id_input}/${j}/g" run0.sh > run1.sh
#sed "s/{outcome_input}/CVD/g" run1.sh > run.sh
#qsub run.sh

done
