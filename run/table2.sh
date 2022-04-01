#!/bin/bash
	
ft="tree"
forms=('-TDPLP' '-TDPBLP' '-TRILP' '-TCFLP')
pcuts=0
prval="0"
table=table2
if test -f "./out/$table.txt"; then
    rm -f "./out/$table.txt"
fi

if test -f "./out/$table.csv"; then
    rm -f "./out/$table.csv"
fi

# Table 2 - LP Relaxation on Trees for TDP, TDP-Benders, TRI, TCF

echo "Instance Formulation Graph_type n m r lp-relaxation runtime num_cuts separation_time" >> ./out/$table.txt

for form in "${forms[@]}";
do 
  for i in {0..74}
  do
    OUTPUT=$(./bin/main -e -u ./dat/${ft}/${ft}-${i}.dat ${form} ${pcuts} ${prval} | tail -1)
    LINE="${ft}-${i}.dat ${form} ${OUTPUT}"
    echo "Solving ${ft}-${i}.dat with ${form}"
    echo ${LINE} >> ./out/$table.txt
  done
done

python ./run/table-constructor.py 2 ./out/$table.txt ./out/$table.csv