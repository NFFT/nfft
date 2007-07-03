hosts=""

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do
hosts="$hosts sfspc$i"
done

for host in $hosts
do
( ssh $host $1 ) &
sleep 1
done
sleep 20
killall ssh
