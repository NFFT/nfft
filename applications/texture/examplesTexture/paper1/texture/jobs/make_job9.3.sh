h=h73
r=r9791
omega=omega9
N=$1
name=omega9.3

jobName=todo.N_$N.$h.$r.$name

echo $h > $jobName
echo $r >> $jobName
echo $omega >> $jobName
echo $N >> $jobName
echo "2 stop_$jobName" >> $jobName
echo >> $jobName
echo $name >> $jobName
