
CC := gcc
CFLAGS :=
OPT := -O3

SRCFILES := gilbert.c

all: gilbert

gilbert: gilbert.c
	$(CC) gilbert.c -o gilbert $(CFLAGS) $(OPT)


