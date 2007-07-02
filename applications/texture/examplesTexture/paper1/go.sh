while [ "`date +%m%d%H`" \< "$1" ]
do
echo Not yet.
date
sleep 60
done

echo Starting.

while [ "`date +%m%d%H`" \< "$2" ]
do
date
./cleanup.sh
./run.sh texture
./run.sh reconstruction
echo sleeping ...
sleep 600
done

echo Cleaning.
./kill.sh
