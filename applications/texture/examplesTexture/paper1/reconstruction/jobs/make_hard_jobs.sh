N=( 5 10 20 40 80 )
N1=( 11 23 41 92 164 308 )
N2=( 128 510 1652 8556 26772 94900 )

# 1 run
for i in 0 1 2 3 4
do
for j in 3 4 5
do
for k in 0 1 2 3 4 5 6 7 8 9
do
./make_job.sh ${N[$i]} ${N1[$j]} ${N2[$j]} 0 $k 1
./make_job.sh ${N[$i]} ${N1[$j]} ${N2[$j]} 1 $k 1
done
done
done
