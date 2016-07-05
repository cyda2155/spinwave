#! /bin/bash
make clean
make
qsub *.pbs
