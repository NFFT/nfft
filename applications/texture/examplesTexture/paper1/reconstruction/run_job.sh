jobdir=jobs
jobName=$1
lockDir=../workinghosts

cat $jobdir/running.$jobName | nice ./main >$jobdir/log.$jobName 2>&1
success=$?

rm $lockDir/`hostname`

# Ready
cd $jobdir
if [ $success -eq 0 ]
then
mv running.$jobName done.${jobName}
else
echo Job failed!
mv running.$jobName failed.$jobName
exit 1
fi


