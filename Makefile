## Makefile 

CC	       = gcc
CC_OPTIONS = -Wall -Wextra -g -v -O3 -std=gnu11 -DNDEBUG
INCLUDES   =  
CFLAGS     =  $(CC_OPTIONS) $(INCLUDES)
LIBS       = -lm

PRIMER_LG_EXE = build/mpd
PRIMER_LG_SRC = src/mem.c src/mpd_lessGreedy.c
PRIMER_LG_OBJ = $(PRIMER_LG_SRC:.c=.o)

PRIMER_MG_EXE  = build/mpd_greedy
PRIMER_MG_SRC  = src/mem.c src/mpd_moreGreedy.c
PRIMER_MG_OBJ  = $(PRIMER_MG_SRC:.c=.o)

POOL_EXE = build/pool_check
POOL_SRC = src/pool_check.c
POOL_OBJ = $(POOL_SRC:.c=.o)

INDEX_EXE = build/index_genome
INDEX_SRC = src/index_genome.c
INDEX_OBJ = $(INDEX_SRC:.c=.o)

PROGS = $(PRIMER_LG_EXE) $(PRIMER_MG_EXE) $(POOL_EXE) $(INDEX_EXE)

all: introduce $(PROGS)
	@echo done.

$(PRIMER_LG_EXE): $(PRIMER_LG_OBJ)
	$(CC) -o $@ $(CFLAGS) $(PRIMER_LG_OBJ) $(LIBS)

$(PRIMER_MG_EXE): $(PRIMER_MG_OBJ)
	$(CC) -o $@ $(CFLAGS) $(PRIMER_MG_OBJ) $(LIBS)

$(POOL_EXE): $(POOL_OBJ)
	$(CC) -o $@ $(CFLAGS) $(POOL_OBJ) $(LIBS)

$(INDEX_EXE): $(INDEX_OBJ)
	$(CC) -o $@ $(CFLAGS) $(INDEX_OBJ) $(LIBS)

introduce:
	@echo "Building..."

clean:
	rm -f *.o

distclean:
	rm -f $(INDEX_EXE) $(PRIMER_LG_EXE) $(PRIMER_MG_EXE) *.o

## end of Makefile
# DO NOT DELETE THIS LINE -- make depend depends on it.
