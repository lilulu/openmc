#!/bin/bash
oldname=20-2G-2fb-expand-2-r000
newname=20-2G-2fb-expand-2

for i in `seq 123 300`;
do
    mkdir $newname-r$i
    cp $oldname/*.xml $newname-r$i/
    cat $newname-r$i/settings.xml | sed s/@@SEED@@/$i/g > $newname-r$i/settings_new.xml
    mv $newname-r$i/settings_new.xml $newname-r$i/settings.xml
done 


