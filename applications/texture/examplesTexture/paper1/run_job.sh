# Host idle?
lockdir=workinghosts
if [ -e $lockdir/`hostname` ]
then
echo `hostname` is already working.
exit 1
else
touch $lockdir/`hostname`
fi

# Run job.
cd $1
./get_job.sh
