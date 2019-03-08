#!/bin/bash
cd /data/users/robertbt/WELibsmolData
for n in {2..20};
do
cd $1${n}
rm $1${n}.log
cd ..
done
