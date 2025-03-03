#!/bin/bash
aux=0
for i in 10 20 30 40 50 60; do
  aux=$i
  aux=$((aux + 10))
  sed -i "4s/$i/$aux/" algwolff.f90
  gfortran -O3 algwolff.f90 -o a1
  ./a1
  if [ ! -d "l$aux" ]; then
    mkdir "l$aux"
  fi
  mv -t "l$aux" fort.* --force
done
sed -i "4s/70/10/" algwolff.f90
