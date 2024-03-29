Note: We use for experiments a dedicated machine using Ubuntu 18.04.4 LTS, an Intel(R) Xeon(R) W-2145 CPU @ 3.70GHz with 16 cores and 512GB RAM, but only a single core is used. Memory measurements are taken through using "/usr/bin/time -v" and reading the "Maximum resident set size" value. The code is written in C++11 and compiled using GCC 7.5.0 with the compiler flags "-O3" and "-std=c++11".

## 1.) PREPARE THE SYSTEM

- We require some modern unix-based system. We tested the following steps on Mac OS X 10.13.6 and Ubuntu 18.04.4.

## 2.) COMPILE THE CODE

- Execute "make joinsampling" from the main directory (where the "Makefile"-file is located)
This will also create the subfolders for the data.

## 3.) DOWNLOAD/PREPARE THE DATA

- TPC-H: The TPC-H benchmark can be downloaded from http://tpc.org and has to be run to generate the csv files.
- DBLP: How to obtain the DBLP dataset is explained on https://anonymous.4open.science/r/cnd-D375. The raw data can be found on http://arnetminer.org/citation.
- TWITTER: The twitter dataset can be obtained from http://an.kaist.ac.kr/traces/WWW2010.html

Each "database" and its schema and views of the queries used in the experiments are defined in a database file using JSON (tpch.json,dblp.json and twitter.json). Currently only csv-type files are supported as data sources, but ODBC support is essentially already in the code.

The ".json" files in the main folder reveal how the data should be named and where everything should be placed:

- In the folder "data/dblp": citation.csv, paper.csv, author.csv, authored.csv, venue.csv
- In the folder "data/tpch/1X": lineitem.tbl, customer.tbl, nation.tbl, orders.tbl,  supplier.tbl, etc
- In the folder "data/tpch/10X": lineitem.tbl, customer.tbl, nation.tbl, orders.tbl,  supplier.tbl, etc
- In the folder "data/tpch/100X": lineitem.tbl, customer.tbl, nation.tbl, orders.tbl,  supplier.tbl, etc
- In the folder "data/twitter": 1_50x.txt, 2_50x.txt, celebrities_profiles.txt

## 4.) EXECUTE THE CODE

- The code receives as a parameter a valid MySQL query except the SQL-like keyword WEIGHTED BY. For instance, the following command (executed from same Folder as "Makefile" file) will execute the TPC-H query WQY with the stream sampler:
```
bin/joinsampling "SELECT * from QY WEIGHTED BY ((e1*(1-d1))*t1*(e2*(1-d2))*t2) LIMIT 1000000 /* db='tpch.json', seed='test123', scalefactor=1 */"
```
- "QY" is here the desired view from the ".json"-database file
- "WEIGHTED BY" is an SQL-like expression followed by the weighting function
- "LIMIT" determines the desired sample size
- "SEED" picks a seed for the sampling process
- "scalefactor=1" chooses the desired scale factor 1 (alternatively 10/100) for TPC-H datasets (ignored for other datasets)

Alternative parameters in the SQL comments:

- "fk=10" will activate the foreign-key economic sampler with a 10x larger uniform sample
- "hash=1" will activate the acyclic economic sampler
- "simplify=1" will activate the cyclic economic sampler
- "inversion=1" will active the inversion sampler
- "ks=1" chooses to perform a K-S-test on the sample
- "write=output.txt" chooses to write the sample to "output.txt"
