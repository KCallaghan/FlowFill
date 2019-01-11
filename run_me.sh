for x in  0.01 
do
echo $x
date > Argentina_startendtimes_$x.txt
filename="Argentina_new_threshold_output_"
mpirun -np 4 ./test $x Argentina_120m.nc $filename$x
date >> Argentina_startendtimes_$x.txt
done

