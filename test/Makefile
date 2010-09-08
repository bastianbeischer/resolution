test: test.cc
	g++ -o test.o -I$(shell root-config --incdir) -I../RES_Event -c test.cc
	g++ -o test $(shell root-config --libs) -L../RES_Event -l RES_Event test.o

clean:
	@echo "Cleaning up test ..."
	@rm -f test.o test