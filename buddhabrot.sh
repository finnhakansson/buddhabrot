#!/bin/sh


MAXIT=1000
XPIX=150
YPIX=100

echo "time nice -n 20 ./buddhabrot --maxit=$MAXIT --xpix=$XPIX --ypix=$YPIX"
time nice -n 20 ./buddhabrot --maxit=$MAXIT --xpix=$XPIX --ypix=$YPIX


