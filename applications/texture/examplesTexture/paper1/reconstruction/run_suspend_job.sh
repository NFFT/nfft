./run_job.sh $1 &

while [ $? -eq 0 ]
do
	sleep 10
	if ./unused.sh
	then 
		killall -18 main
	else
		killall -19 main
	fi
done
