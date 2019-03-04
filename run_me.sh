for x in  0.2
do
echo $x
date > temp_startendtimes_$x.txt
filename="Sangamon_text_"
mpirun -np 4 ./test $x Sangamon_15_m.nc $filename$x
date >> temp_startendtimes_$x.txt
done

#0.001 0.01 0.05 0.1 0.2