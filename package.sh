#!/bin/bash

make -j 4

mkdir bin

mv bgreat bin

mv numbersToSequences bin

tar -czvf bin.tar.gz bin;


echo The end !;

