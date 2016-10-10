
.PHONY: all clean

all: poisson-rect.out

%.out: %.o
	$(CC) -o $@ $<

%.o: %.c
	$(CC) -c -o $@ $<

clean:
	git clean -fX .

