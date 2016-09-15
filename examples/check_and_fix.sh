#!/bin/bash
#name=2dbeavrs-full-2/2fb-expand
name=2dbeavrs-full-1-4m/2fbc-rd-f20
set -e
#for i in 100-2
#for i in 108 133
for i in `seq 100 104`;
do
    #if ! test -f "$name/r$i/statepoint.500.binary"; then 
    if grep -q CANCELLED $name/r$i/*.err; then
    #if true; then
        echo found case $name/r$i killed
        if test -f "$name/r$i/fs.dat"; then
            mv $name/r$i/fs.dat $name/r$i/fs.old.dat
        fi
        if test -f "$name/r$i/log.txt"; then
            mv $name/r$i/log.txt $name/r$i/log.old.txt
        fi
        if [ `ls -1 $name/r$i/slurm.*.out 2>/dev/null | wc -l ` -gt 0 ]; then
            for f in $name/r$i/slurm.*.out; do mv $f $f.old; done
            for f in $name/r$i/slurm.*.err; do mv $f $f.old; done
        fi
        (cd "$name/r$i" && pwd && sbatch run.slurm)
    fi
done 


