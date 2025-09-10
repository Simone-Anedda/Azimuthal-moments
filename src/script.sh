#!/bin/bash

FILE=$1
END="$2"
OUTFILE="./output/res-ref-fit-2018-bT-space-16-10-23-MSHT20nlo_as118-DEHSS14NLO-CT3evo-bootstrap-replicas.dat"
#echo $END

function escape_slashes {     sed 's/\//\\\//g' ; }

function change_line {
     local OLD_LINE_PATTERN=$(grep 'index' $1);
#     echo $OLD_LINE_PATTERN;
#     local seed=$2
#      echo "seed="$2
     local NEW_LINE="index "$2;
#     local NEW_LINE="set seed = "$(( $RANDOM % 100 + 1 ));
#     local FILE=$1;
#     local (NEW=$NEW_LINE | escape_slashes);
     local NEW=$(echo "${NEW_LINE}" | escape_slashes)
#      echo "new="$NEW; 
#      echo $FILE;
     sed -i '/'"${OLD_LINE_PATTERN}"'/s/.*/'"${NEW}"'/' "${FILE}";
#      grep "set seed" $FILE
}

#change_line $FILE

#./ho_cluster < $FILE &
i=0
touch $OUTFILE
# time ./exe/fit_main $(cat $FILE) &
# echo $(wc -l < $OUTFILE)
# length=$(cat "$OUTFILE" | wc -l)
# for ((i = 0 ; i <= $END ; i++)); do
while (( $(cat "$OUTFILE" | wc -l) < $END )); do
#   echo $i
#  then
    time ./exe/fit_main $(cat $FILE) &
    ((i+=1))
#     length=$(cat "$OUTFILE" | wc -l)
#     echo $length $i
    change_line $FILE $i
    sleep 2;
#  fi
done
