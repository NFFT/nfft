N=( 05 10 20 40 80 )
N1=( 11 23 41 92 164 308 )
N2=( 128 510 1652 8556 26772 94900 )

# 10 runs

for i in 0 1 2 3 4
do
./make_job.sh ${N[$i]} ${N1[0]} ${N2[0]} 0 0 10
./make_job.sh ${N[$i]} ${N1[0]} ${N2[0]} 1 0 10
done

for i in 0 1 2 3
do
./make_job.sh ${N[$i]} ${N1[1]} ${N2[1]} 0 0 10
./make_job.sh ${N[$i]} ${N1[1]} ${N2[1]} 1 0 10
done

# 5 runs

./make_job.sh ${N[4]} ${N1[1]} ${N2[1]} 0 0 5
./make_job.sh ${N[4]} ${N1[1]} ${N2[1]} 0 5 5
./make_job.sh ${N[4]} ${N1[1]} ${N2[1]} 1 0 5
./make_job.sh ${N[4]} ${N1[1]} ${N2[1]} 1 5 5

for i in 0 1 2 3 4
do
./make_job.sh ${N[$i]} ${N1[2]} ${N2[2]} 0 0 5
./make_job.sh ${N[$i]} ${N1[2]} ${N2[2]} 0 5 5
./make_job.sh ${N[$i]} ${N1[2]} ${N2[2]} 1 0 5
./make_job.sh ${N[$i]} ${N1[2]} ${N2[2]} 1 5 5
done
