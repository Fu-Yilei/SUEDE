#!/bin/sh

pip install ncbi-genome-download, progressbar
conda install ete3, pandas, seaborn, blast
git clone https://github.com/marbl/parsnp.git 
cd parsnp 
./build_parsnp_linux.sh
