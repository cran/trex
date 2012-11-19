#!/bin/bash


## write table11 to a file, then run trex with threshold=2 on it
echo "2 46" >> table11.txt
echo "1 46" >> table11.txt
echo "8 98" >> table11.txt

trex -table table11.txt -threshold 2 > test11.out


## write table10 to a file, then run trex with threshold=2 on it
echo "8 92" >> table10.txt
echo "0  0" >> table10.txt
echo "2 98" >> table10.txt

trex -table table10.txt -threshold 2 > test10.out


## run trex with threshold=4 on table10
trex -table table10.txt -threshold 4 > test10.4.out