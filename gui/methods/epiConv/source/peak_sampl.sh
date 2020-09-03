#!/bin/bash
#usage: peak_sampl.sh <peakfile> <sample times> <fraction of sampled peaks> <random seed>

NPEAK=$(wc -l $1|gawk '{print $1}')
NCOL=$[$2]
echo "$NPEAK $NCOL"
gawk 'BEGIN{ORS="";srand('"$4"')}{for(i=1;i<='"$2"';i++){if(rand()<'"$3"')print $4,i"\n"}}' $1
