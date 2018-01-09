#!/bin/sh
tar vxf samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make
cd ..
make
gcc -O2 chrLenExtract.c -o chrLenExtract
chmod +x ./chrLenExtract
chmod +x ./BS-Snper
