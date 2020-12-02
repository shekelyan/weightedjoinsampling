1.) PREPARE THE DATASETS

TPC-H: The TPC-H benchmark can be downloaded from http://tpc.org and has to be run to generate the csv files.
DBLP: How to obtain the DBLP dataset is explained by one of our authors on https://github.com/qingzma/cnd. 
TWITTER: The twitter dataset can be obtained from http://arnetminer.org/citation.

2.) MODIFY THE DATABASE FILES

The ".json" files in the main folder reveal how the data should be named and where everything should be placed.

3.) COMPILE THE CODE

Execute "make joinsampling" from the directory where the "Makefile"-file is located.

4.) EXECUTE THE CODE

The following command will execute the TPC-H query WQY with the stream sampler.:

bin/joinsampling "SELECT * from QY WEIGHTED BY ((e1*(1-d1))*t1*(e2*(1-d2))*t2) LIMIT 1000000 /* db='tpch.json', seed='test123', scalefactor=1 */"

To activate the foreign-key economic sampler one adds the parameter "fk=10" in the SQL comment. The 10 signifies how much larger of a sample it collects. To activate the many-to-many economic sampler one adds the parameter "hash=1" in the SQL comment. To activate the cyclic economic sampler one adds the parameter "simplify=1" in the SQL comment. To activate the index-based sampler one uses adds parameter "index=1" in the SQL comment. To activate the inversion sampler one uses the parameter "inversion=1". For KS-testing one adds the parameter "ks=1".
