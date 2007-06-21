jobdir=jobs
jobName=$1
lockDir=../../workinghosts

cat $jobdir/running.$jobName | nice ./main >$jobdir/log.$jobName 2>&1
success=$?

# Ready
cd $jobdir
if [ $success -eq 0 ]
then

rm $lockDir/`hostname`
mv running.$jobName done.$jobName
rm log.$jobName

else

rm $lockDir/`hostname`
echo Job failed!
mv running.$jobName failed.$jobName
exit 1

fi
