#!/bin/bash
cd /dfs3/pub/robertbt/WELibsmolData/

for n in {21..100}
do
for runNo in {1..3}
do
cd ${1}${n}runNo${runNo}
if [ -e savestate.txt ]
then
qsub ./${1}${n}runNo${runNo}.pub
fi
cd ..
done
done