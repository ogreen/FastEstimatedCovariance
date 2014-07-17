INC = -I /usr/local/include

CFLAGS= -g -msse -mfused-madd 
CFLAGS+=-Iinclude -fopenmp

LDLIBS=-lm -lrt

SRC = Main.c BuildCovMat.c BuildCombinationsTable.c RunCombination.c timing_util.c

SRCALL= $(SRC)

all: $(SRCALL)
	gcc -std=gnu9x  $(CFLAGS) -O2 -o $@ $^ $(LDFLAGS) $(LDLIBS)

.c.o:
	gcc -fopenmp -O0 -std=gnu9x -Wall $(INC) -c $<

clean:
	rm -f all $(OBJS)
