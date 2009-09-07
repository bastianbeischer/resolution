LIB = libBlobel.so

FC = gfortran
LD = gfortran

FFLAGS += -fPIC

all: $(LIB)

clean:
	@echo "Cleaning up blobel ..."
	@rm -f *.o $(LIB)

$(LIB): blobel.o
	$(LD) -o $(LIB) blobel.o -shared
	@echo "DONE."
