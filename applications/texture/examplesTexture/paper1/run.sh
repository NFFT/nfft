workDir=trunk/applications/texture/examplesTexture/paper1
workScript=run_job.sh
../texture_util/forallhosts.sh "cd $workDir && ./$workScript $1"
