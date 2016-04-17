#!/bin/bash
imodel=1
#01 : cuprate
#02 : cuprate2l
case $imodel in
  1)
  dir='cuprate'
  u=6.0 
  vnn=-1.6
  t=-0.00
  mu=-1.0
  la=0.01
  nlform=9
  ng=10
  ;;
  2)
  dir='cuprate2l'
  u=2.0
  vnn=0.0
  mu=9.26
  ;;
esac
mkdir -p $dir
cd $dir
bcsflow='f'
dataready='f'
quicksearch='f'
skipinterorbitpair='f'
appendbreakpoint='f'
appendkmesh='f'
useprojectors='f'
appendprojectors='f'
periodiclayer='f'
nw=400
wmax=100
wmin=1.e-3
wir=1.e-4
diverge=80
vl2l=0
vbonding=0
vantibonding=0
jnn=0
dir1='u'$u'v'$vnn'r'$la't'$t'mu'$mu'nlform'$nlform'ng'$ng
mkdir -p $dir1
cd $dir1
cp ../../smfrg_lee.out .
echo $imodel > imodel.input
echo $bcsflow > bcsflow.input
echo $dataready > dataready.input
echo $quicksearch > quicksearch.input
echo $appendbreakpoint > appendbreakpoint.input
echo $appendkmesh > appendkmesh.input
echo $useprojectors > useprojectors.input
echo $appendprojectors > appendprojectors.input
echo $nw > nw.input
echo $wmax > wmax.input
echo $wmin > wmin.input
echo $wir > wir.input
echo $diverge > diverge.input
echo $u > u.input
echo $vnn > vnn.input
echo $mu > mu.input
echo $vl2l > vl2l.input
echo $vbonding > vbonding.input
echo $vantibonding > vantibonding.input
echo $jnn > jnn.input
echo $la > la.input
echo $t > t.input
echo $nlform > nlform.input
echo $ng > ng.input
# echo ' ' > output
