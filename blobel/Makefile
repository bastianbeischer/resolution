LIB = libBlobel.so

FC = gfortran
LD = gfortran

FFLAGS += -fPIC

all: $(LIB)

clean:
		@rm -f *.o $(LIB)

$(LIB): blobel.o
	@echo "Creating shared lib libBlobel.so"
	$(LD) -o $(LIB) blobel.o -shared
	@echo "DONE."
	@echo
