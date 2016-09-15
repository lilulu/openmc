#!/bin/bash
name=20-2G-2fbc-rd-expand-2

for i in `seq 200 299`;
do
    mkdir $name/r$i
    cp $name/r000/*.xml $name/r$i/
    cat $name/r$i/settings.xml | sed s/@@SEED@@/$i/g > $name/r$i/settings_new.xml
    mv $name/r$i/settings_new.xml $name/r$i/settings.xml
done 


