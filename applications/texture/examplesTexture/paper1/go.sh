while [ "`date +%d%H%M%S`" \< "$1" ]
do
echo Not yet.
date
sleep 60
done

echo Starting.

while [ "`date +%d%H%M%S`" \< "$2" ]
do
date
./cleanup.sh
./run.sh reconstruction
echo sleeping ...
sleep 300
done

echo Cleaning.
./kill.sh
