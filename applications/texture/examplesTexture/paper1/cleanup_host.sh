lockdir=workinghosts
jobdir=reconstruction/jobs

procs=`ps x`
mainprocs=`echo $procs | grep ./main`
echo $mainprocs

if [ "$mainprocs" = "" ]
then
if [ -e "$lockdir/`hostname`" ]
then
echo `hostname` has been aborted.

rename running todo $jobdir/running.*.`hostname`
rm -v "$lockdir/`hostname`"

fi
fi
