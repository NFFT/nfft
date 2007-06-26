while [ 0 ]
do
date
./cleanup.sh
./run.sh reconstruction
echo sleeping ...
sleep 1000
done
