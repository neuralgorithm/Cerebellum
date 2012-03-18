CC = gcc
CFLAGS = -O3 -Wall

MT = mt19937ar-cok

all: okr pfinit si ri wdist

clean:
	rm -f *.o okr pfinit si ri wdist

okr: okr.o $(MT).o
	$(CC) $(CFLAGS) -o $@ $@.o $(MT).o -lm -lgsl -lblas

pfinit: pfinit.o
	$(CC) $(CFLAGS) -o $@ $@.o -lm

si: si.o
	$(CC) $(CFLAGS) -o $@ $@.o -lm -lgd

ri: ri.o
	$(CC) $(CFLAGS) -o $@ $@.o -lm -lgd

wdist: wdist.o
	$(CC) $(CFLAGS) -o $@ $@.o -lm
