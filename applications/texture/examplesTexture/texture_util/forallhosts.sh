hosts=""
if [ -z "$2" ]
then
	delay=5
else
	delay=$2
fi
if [ -z "$3" ]
then
	final=20
else
	final=$3
fi

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do
hosts="$hosts schmalzm@sfspc$i"
done

for host in $hosts
do
( ssh $host $1 ) &
sleep $delay
done
sleep $final
killall ssh
