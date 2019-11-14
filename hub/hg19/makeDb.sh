#! /usr/bin/env bash

rm -f trackDb.txt 

bash dms/makeDmsDb.sh > dms/trackDb.txt
bash enzymatic_probing/makeEpDb.sh > enzymatic_probing/trackDb.txt
bash rnastructurome/makeTrackdb.sh > rnastructurome/trackDb.txt
bash editing/makeTrackdb.sh > editing/trackDb.txt
bash icSHAPE/makeTrackdb.sh > icSHAPE/trackDb.txt

dbfiles=$(find . -maxdepth 2 -name 'trackDb.txt')
cat tDb.tmp $dbfiles >> trackDb.txt

