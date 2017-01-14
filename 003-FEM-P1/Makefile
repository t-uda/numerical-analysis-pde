
CFLAGS=-Wall -O3
LDFLAGS=-lm

.PHONY: all clean
.PRECIOUS: %.c %.o

all: fem.out rectmesh.out genmesh.out

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.out: %.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

fem.out: fem.o mesh.o vector.o
rectmesh.out: rectmesh.o
genmesh.out: genmesh.o mesh.o vector.o

clean:
	$(RM) *.o *.out

