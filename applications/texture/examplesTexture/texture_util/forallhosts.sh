max64=48
max21=10
hosts=""
user=schmalzm

i=1
while [ $i -le $max64 ]
do
if [ $i -le 9 ]
then
number=0$i
else
number=$i
fi
hosts="$hosts $user@64pc$number"
let i=$i+1
done

hosts="$hosts $user@64pc66"

i=1
while [ $i -le $max21 ]
do
if [ $i -le 9 ]
then
number=0$i
else
number=$i
fi
hosts="$hosts $user@21pc$number"
let i=$i+1
done

for host in $hosts
do
ssh $host $1
done
