LIB = libMillepede.so

FC = gfortran
LD = gfortran

FFLAGS += -fPIC

all: $(LIB)

clean:
	@echo "Cleaning up millepede ..."
	@rm -f *.o $(LIB)

$(LIB): millepede1.o
	$(LD) -o $(LIB) millepede1.o -shared
