#!/bin/bash
#name=2dbeavrs-full-1-4m/2fbc-rd-p
name=2dbeavrs-full-1-4m/2fbc-rd-p
set -e
for i in `seq 101 104`;
do
    mkdir $name/r$i
    cp $name/r000/run.slurm $name/r$i/run.slurm
    cat $name/r$i/run.slurm | sed s/@@SEED@@/$i/g > $name/r$i/run-new.slurm
    mv $name/r$i/run-new.slurm $name/r$i/run.slurm
    cp $name/r000/*.xml $name/r$i/
    cat $name/r$i/settings.xml | sed s/@@SEED@@/$i/g > $name/r$i/settings_new.xml
    mv $name/r$i/settings_new.xml $name/r$i/settings.xml
    (cd "$name/r$i" && pwd && sbatch run.slurm)
done 


