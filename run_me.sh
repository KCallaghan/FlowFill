input="user_inputs.txt"
counter=0
while IFS= read -r var
do
  tmp=${var#*: }
  counter=$((counter+1))
  if [ $counter == 1 ]
  then
    processors="$tmp"
  elif [ $counter == 2 ]
  then
    runfile="$tmp"
  elif [ $counter == 3 ]
  then
    topofile="$tmp"
  elif [ $counter == 4 ]
  then
    cols="$tmp"
  elif [ $counter == 5 ]
  then
    rows="$tmp"
  elif [ $counter == 6 ]
  then
    threshold="$tmp"
  elif [ $counter == 7 ]
  then
    outfile="$tmp"
  elif [ $counter == 8 ]
  then
    runoff_bool="$tmp"
  elif [ $counter == 9 ]
  then
    runoff="$tmp"
  elif [ $counter == 10 ]
  then
    runoff_file="$tmp"
  elif [ $counter == 11 ]
  then
    ties="$tmp"
  fi
done < "$input"


date > temp_check_startendtimes_$x.txt
mpirun -np $processors ./$runfile $runoff $topofile $rows $cols $threshold $outfile $runoff_bool $runoff_file $ties
date >> temp_check_startendtimes_$x.txt


#0.001 0.01 0.05 0.1 0.2
