#!/bin/bash

ft="series-parallel"
forms=('-FLOWLP')
pcuts_list=(0 100)
prval="100.0"

table=table3
if test -f "./out/$table.txt"; then
    rm -f "./out/$table.txt"
fi

# Table 3 - Effect of DP cuts on the LP bound of the flow formulation

echo "Instance Formulation Graph_type n m r lp-relaxation runtime cut_percentage num_cuts" >> ./out/$table.txt

for form in "${forms[@]}";
do
  for pcuts in ${pcuts_list[@]};
  do
    for i in {0..2} {5..7} {10..12} {16..18} {24..26} {33..35};
    do
      if [[ $pcuts -gt 0 ]]
      then
        SFX="${form}+DP "
      else
        SFX="${form} "
      fi
      OUTPUT=$(./bin/main -e -u ./dat/${ft}/${ft}-${i}.dat ${form} ${pcuts} ${prval} | tail -1)
      LINE="${ft}-${i}.dat ${SFX}${OUTPUT}"
      echo "Solving ${ft}-${i}.dat with ${form}"
      echo ${LINE} >> ./out/$table.txt
    done
  done
done

