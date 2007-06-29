N=( 05 10 20 40 80 )
N1=( 011 023 041 092 164 308 )
N2=( 00034 00156 00460 02248 06974 23898 )

# 10 runs

for i in 0 1 2 3 4
do
./make_job2.sh ${N[$i]} ${N1[0]} ${N2[0]} 0 0 10
./make_job2.sh ${N[$i]} ${N1[0]} ${N2[0]} 1 0 10
done

for i in 0 1 2 3
do
./make_job2.sh ${N[$i]} ${N1[1]} ${N2[1]} 0 0 10
./make_job2.sh ${N[$i]} ${N1[1]} ${N2[1]} 1 0 10
done

# 5 runs

./make_job2.sh ${N[4]} ${N1[1]} ${N2[1]} 0 0 5
./make_job2.sh ${N[4]} ${N1[1]} ${N2[1]} 0 5 5
./make_job2.sh ${N[4]} ${N1[1]} ${N2[1]} 1 0 5
./make_job2.sh ${N[4]} ${N1[1]} ${N2[1]} 1 5 5

for i in 0 1 2 3 4
do
./make_job2.sh ${N[$i]} ${N1[2]} ${N2[2]} 0 0 5
./make_job2.sh ${N[$i]} ${N1[2]} ${N2[2]} 0 5 5
./make_job2.sh ${N[$i]} ${N1[2]} ${N2[2]} 1 0 5
./make_job2.sh ${N[$i]} ${N1[2]} ${N2[2]} 1 5 5
done

# 1 run
for i in 0 1 2 3 4
do
for j in 3 4 5
do
for k in 0 1 2 3 4 5 6 7 8 9
do
./make_job2.sh ${N[$i]} ${N1[$j]} ${N2[$j]} 0 $k 1
./make_job2.sh ${N[$i]} ${N1[$j]} ${N2[$j]} 1 $k 1
done
done
done
