CC = g++
FLAGS = -O2 
LDFLAGS = -O2
LIBS = -m64 -I./samtools-0.1.19/ -L./samtools-0.1.19/ -lbam -lz -lpthread

TARGET = bsSnp
TOCOMPILE = main.o sam_funcs.o hash_funcs.o chrome_funcs.o 

all: ${TOCOMPILE}
	${CC} $(LDFLAGS) -o $(TARGET) ${TOCOMPILE} ${LIBS}

.c.o:
	$(CC) ${FLAGS} ${LIBS} -c $*.c

clean:
	rm -f *.o $(TARGET)

