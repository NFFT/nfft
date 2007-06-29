NList="000 002 004 008 016 032 064 128"

for N in $NList
do
./make_job.sh $N
done
