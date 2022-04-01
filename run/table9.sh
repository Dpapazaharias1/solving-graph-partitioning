#!/bin/bash

fts=('series-parallel' 'toroidal-grid' 'planar-grid')
forms=('-TRI' '-FLOW' '-PATH' '-FLOW+')

table=table9
if test -f "./out/$table.txt"; then
    rm -f "./out/$table.txt"
fi
# Table 9: Formulation strength of IP formulations on special graph structures

#echo "Instance -TRI graph_name n m r best_obj mip_gap runtime b&b nodes" >> ./out/$table.txt
#echo "Instance -FLOW graph_name n m r best_obj mip_gap runtime b&b_nodes cut_type prob_cut r_pct_trees_cut num_cuts separation_time" >> ./out/$table.txt
#echo "Instance -PATH graph_name n m r best_obj mip_gap runtime b&b_nodes num_lazy num_user lazy_time user_time -tree_cut_type(not used)" >> ./out/$table.txt
#echo "Instance -FLOW+ graph_name n m r best_obj mip_gap runtime b&b_nodes cut_type prob_cut r_pct_trees_cut num_cuts separation_time" >> ./out/$table.txt

echo "Instance Formulation Graph_type n m r best_obj mip_gap runtime b&b_nodes cut_type prob_cut r_pct_trees_cut num_cuts separation_time x x" >> ./out/$table.txt

for ft in "${fts[@]}";
do
    if [ $ft == 'series-parallel' ]
    then 
        for i in {0..42};
        do
            for form in "${forms[@]}";
            do
                if [ $form == '-FLOW+' ]
                then
                    OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat ${form} 50 1.00 -prim | tail -1)
                else
                    OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat ${form} 0 0 | tail -1)
                fi
                LINE="${ft}-${i}.dat ${form} ${OUTPUT}"
                echo "Solving ${ft}-${i}.dat with ${form}"
                echo ${LINE} >> ./out/$table.txt
            done
        done
    fi

    if [ $ft == 'toroidal-grid' ]
    then
        for i in {0..52};
        do
            for form in "${forms[@]}";
            do
                if [ $form == '-FLOW+' ]
                then
                    OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat ${form} 50 1.00 -prim | tail -1)
                else
                    OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat ${form} 0 0 | tail -1)
                fi
                LINE="${ft}-${i}.dat ${form} ${OUTPUT}"
                echo "Solving ${ft}-${i}.dat with ${form}"
                echo ${LINE} >> ./out/$table.txt
            done
        done
    fi

    if [ $ft == 'planar-grid' ]
    then 
    for i in {0..111};
    do
        for form in "${forms[@]}";
        do
            if [ $form == '-FLOW+' ]
            then
                OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat ${form} 50 1.00 -prim | tail -1)
            else
                OUTPUT=$(./bin/main -e -w ./dat/${ft}/${ft}-${i}.dat ${form} 0 0 | tail -1)
            fi
            LINE="${ft}-${i}.dat ${form} ${OUTPUT}"
            echo "Solving ${ft}-${i}.dat with ${form}"
            echo ${LINE} >> ./out/$table.txt
        done
    done
    fi
    
done


