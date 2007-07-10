jobdir=jobs
lockdir=../../workinghosts

# Choose job.
cd $jobdir

for count in 1 2 3 4 5
do
jobNameOld=`ls -1 todo.* 2>/dev/null | head -1`
if ! ls todo.* 2>/dev/null 1>/dev/null
then
rm $lockdir/`hostname`
echo No jobs left!
exit 1
fi


jobNameOld=${jobNameOld#todo.}

# Block job.
jobName=$jobNameOld.`hostname`
mv todo.$jobNameOld running.$jobName

run=$?
if [ $run -eq 0 ]
then
break
fi

done

if [ $run -ne 0 ]
then
rm $lockdir/`hostname`
echo Unable to run job!
exit 1
fi

cd ..

# Run job.
echo Running job $jobName.
screen -d -m ./run_suspend_job.sh $jobName
