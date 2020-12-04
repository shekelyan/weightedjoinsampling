joinsampling:
	mkdir -p bin
	mkdir -p data/tpch/1x
	mkdir -p data/tpch/10x
	mkdir -p data/tpch/100x
	mkdir -p data/dblp
	mkdir -p data/twitter
	$(CXX) -w -g -Isrc -O3 -std=c++11 -o bin/joinsampling src/joinsampling.cpp

.PHONY:
	clean
clean:
	rm -f -r bin/*.dSYM

