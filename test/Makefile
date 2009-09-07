test: test.cc
	g++ -o test.o -I$(ROOTSYS)/include -I../RES_Event -c test.cc
	g++ -o test $(shell root-config --libs) -L../RES_Event -l RES_Event test.o

clean:
	rm -f test.o test