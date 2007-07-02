lockdir=workinghosts
jobdir1=reconstruction/jobs
jobdir2=texture/jobs
host=`hostname`

procs=`ps x`
mainprocs1=`echo $procs | grep ./main`
mainprocs2=`echo $procs | grep ./test`
jobs1=`ls "${jobdir1}/running.*.${host}" 2>/dev/null`
jobs2=`ls "${jobdir2}/running.*.${host}" 2>/dev/null`

echo `hostname`

if [ "$mainprocs1" = "" ] && [ "$jobs1" != "" ]
then
echo `hostname` has been aborted.

~/bin/rename running todo $jobdir1/running.*.`hostname`

fi

if [ "$mainprocs2" = "" ] && [ "$jobs2" != "" ]
then
echo `hostname` has been aborted.

~/bin/rename running todo $jobdir2/running.*.`hostname`

fi

if [ "$mainprocs1" = "" ] && [ "$mainprocs2" = "" ] && [ -e workinghosts/`hostname` ]
then
rm -v workinghosts/`hostname`
fi
