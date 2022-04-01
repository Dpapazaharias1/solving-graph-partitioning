#!/bin/bash

ft="series-parallel"
forms=('-TCFIP' '-BLOCKIP')
pcuts_list=100
prval="100.0"

table=table4
if test -f "./out/$table.txt"; then
    rm -f "./out/$table.txt"
fi

# Table 4 - Comparison of TCF IP and Block decomposition IP formulations

echo "Instance Formulation graph_type n m r best_obj mip_gap runtime cut_type num_lazy_cuts lazy_separation_time" >> ./out/$table.txt

for form in "${forms[@]}";
do
  for pcuts in ${pcuts_list[@]};
  do
    for i in {0..2} {5..7} {10..12} {16..18} {24..26} {33..35};
    do
      OUTPUT=$(./bin/main -e -u ./dat/${ft}/${ft}-${i}.dat ${form} ${pcuts} ${prval} | tail -1)
      LINE="${ft}-${i}.dat ${form} ${OUTPUT}"
       echo "Solving ${ft}-${i}.dat with ${form}"
    echo ${LINE} >> ./out/$table.txt
    done
  done
done

python ./run/table-constructor.py 4 ./out/$table.txt ./out/$table.csv


