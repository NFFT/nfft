./run_job.sh $1 &

while [ $? -eq 0 ]
do
	sleep 10
	users=`who | wc -l`
	me=`who | grep schmalzm | wc -l`
	root=`who | grep root | wc -l`
	let foreign=$users-$me-$root
	if [ $foreign -eq 0 ]
	then 
		killall -18 test
	else
		killall -19 test
	fi
done
