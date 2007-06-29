lockdir=workinghosts
jobdir1=reconstruction/jobs
jobdir2=texture/jobs

procs=`ps x`
mainprocs=`echo $procs | grep ./main`
mainprocs=$mainprocs`echo $procs | grep ./test`

echo `hostname`

if [ "$mainprocs" = "" ]
then
if [ -e "$lockdir/`hostname`" ]
then
echo `hostname` has been aborted.

rename running todo $jobdir1/running.*.`hostname`
rename running todo $jobdir2/running.*.`hostname`
rm -v "$lockdir/`hostname`"

fi
fi
