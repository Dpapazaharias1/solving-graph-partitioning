#!/bin/bash

ft="watts-strogatz"
forms=('-FLOWLP')
pcuts_list=(50 100)
prval_list=('1.00' '1.25' '1.50')

table=table5-6
if test -f "./out/$table.txt"; then
    rm -f "./out/$table.txt"
fi

# Table 5-6 - Comparison of heuristically separated prim and knapsack user cuts

echo "Instance Formulation Graph_type n m r best_obj mip_gap runtime b&b_nodes cut_type prob_cut r_pct_trees_cut num_cuts separation_time" >> ./out/$table.txt

for i in {0..80};
do
    OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat -FLOW 0 0 | tail -1)
    LINE="${ft}-${i}.dat -FLOW ${OUTPUT}"
    echo "Solving ${ft}-${i}.dat with -FLOW"
    echo ${LINE} >> ./out/$table.txt
done
#0-80
for i in {0..80};
do
    for pcuts in ${pcuts_list[@]};
    do
        for prval in ${prval_list[@]};
        do
        if [ $prval == "1.00" ]
        then
            OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat -FLOW+ ${pcuts} ${prval} -prim | tail -1)
            LINE="${ft}-${i}.dat -FLOW+ ${OUTPUT}"
            echo "Solving ${ft}-${i}.dat with -FLOW+"
    echo ${LINE} >> ./out/$table.txt
        else
            OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat -FLOW+ ${pcuts} ${prval} -knap | tail -1)
            LINE="${ft}-${i}.dat -FLOW+ ${OUTPUT}"
            echo "Solving ${ft}-${i}.dat with -FLOW+"
            echo ${LINE} >> ./out/$table.txt
        fi
        done
    done
done


python ./run/table-constructor.py 5 ./out/$table.txt ./out/table5.csv
python ./run/table-constructor.py 6 ./out/$table.txt ./out/table6.csv