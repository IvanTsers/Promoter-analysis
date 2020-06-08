#!/bin/bash
#

echo '    Please enter the number of threads:'
read thread
echo '    Thread number is' $thread'.'

mkdir MAST_output

cd ./PlantPAN_meme_motifs

for filename in *.meme; do
 (
  echo
  echo '    Processing file' $filename'...'
  mast $filename ../Promoters.fa -oc ../MAST_output/$filename -mf $filename
 ) &

   # allow only to execute $thread jobs in parallel (one of them as the 'main' job and the rest as 'background' jobs)
   if [[ $(jobs -r -p | wc -l) -gt $(($thread-1)) ]]; then
   # wait only for the first (the 'main') job
      wait -n
   fi
done

# wait for pending jobs
wait
