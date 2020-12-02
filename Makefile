joinsampling:
	$(CXX) -w -g -Isrc -O3 -std=c++11 -o bin/joinsampling src/joinsampling.cpp

.PHONY:
	clean
clean:
	rm -f -r bin/*.dSYM

