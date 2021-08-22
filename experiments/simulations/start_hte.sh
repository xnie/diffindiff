#! /bin/bash

nonconstsetups=('A' 'B' 'C' 'D')
ns=(1000 2000)
ps=(6 12)
rep=200


for ((i1=0; i1<${#nonconstsetups[@]}; i1++))
do
for ((i2=0; i2<${#ns[@]}; i2++))
do
for ((i3=0; i3<${#ps[@]}; i3++))
do
  setup=${nonconstsetups[$i1]}
  n=${ns[$i2]}
  p=${ps[$i3]}

  fnm="logging/progress-$setup-$n-$p-$rep.out"
  echo $fnm

  Rscript run_sims_hte.R $setup $n $p $rep 2>&1 | tee $fnm &
done
done
done
