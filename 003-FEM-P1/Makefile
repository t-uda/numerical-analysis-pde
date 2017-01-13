
CFLAGS=-Wall -O3
LDFLAGS=-lm

.PHONY: all clean
.PRECIOUS: %.c %.o

all: fem.out genmesh.out

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.out: %.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	$(RM) *.o *.out

