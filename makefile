CC = gcc
CFLAGS = $(shell pkg-config --cflags gtk+-3.0)
LIBS = $(shell pkg-config --libs gtk+-3.0)

TARGETS = simplex

all: $(TARGETS)

simplex: simplex.c
	$(CC) simplex.c -o simplex $(CFLAGS) $(LIBS) -lm

clean:
	rm -f simplex
