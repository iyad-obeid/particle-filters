# --- macros
CC=gcc
CFLAGS=  -O3 -I /usr/local/include
OBJECTS = createData.o oMatrixPrimitives.o liuPF.o
LIBS = -lgsl -lgslcblas -lm -L/usr/local/lib
LIBS2 = -lm -L/usr/local/lib


# --- targets
all: liuPF liuVPF createSingleInputFile liuPFx liuVPFx

liuPF: $(OBJECTS)
	$(CC) -o liuPF $(LIBS) $(OBJECTS)

liuVPF: createData.o oMatrixPrimitives.o liuVPF.o oMatrixPrimitives.o
	$(CC) -o liuVPF $(LIBS) createData.o oMatrixPrimitives.o liuVPF.o

liuPFx: liuPFx.o createData.o oMatrixPrimitives.o
	$(CC) -o liuPFx $(LIBS) liuPFx.o createData.o oMatrixPrimitives.o

liuVPFx: createData.o oMatrixPrimitives.o liuVPFx.o oMatrixPrimitives.o
	$(CC) -o liuVPFx $(LIBS) createData.o oMatrixPrimitives.o liuVPFx.o

createSingleInputFile: createSingleInputFile.o createData.o oMatrixPrimitives.o
	$(CC) -o createSingleInputFile $(LIBS) createData.o oMatrixPrimitives.o createSingleInputFile.o

liuVPF.o: liuVPF.c
	$(CC) $(CFLAGS) -c liuVPF.c

liuVPFx.o: liuVPFx.c
	$(CC) $(CFLAGS) -c liuVPFx.c

liuPF.o: liuPF.c
	$(CC) $(CFLAGS) -c liuPF.c

liuPFx.o: liuPFx.c
	$(CC) $(CFLAGS) -c liuPFx.c
        
createSingleInputFile.o: createSingleInputFile.c
	$(CC) $(CFLAGS) -c createSingleInputFile.c

createData.o: createData.c
	$(CC) $(CFLAGS) -c createData.c
       
oMatrixPrimitves.o: oMatrixPrimitives.c
	$(CC) $(CFLAGS) -c oMatrixPrimitives.c

foo: oMatrixPrimitives.o foo.o
	$(CC) -o foo $(LIBS) oMatrixPrimitives.o foo.o

foo.o: foo.c
	$(CC) $(CFLAGS) -c foo.c

# --- remove binary and executable files
clean:
	rm -f liuPF $(OBJECTS)
