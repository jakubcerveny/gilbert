
CC := gcc
RUSTC := rustc
CFLAGS :=
OPT := -O3
#OPT := -g
RUSTOPT := -O

all: gilbert3dpp-c

gilbert3dpp-c: gilbert3dpp.c
	$(CC) gilbert3dpp.c -o gilbert3dpp-c $(CFLAGS) $(OPT)

.PHONY: clean all

clean:
	rm -f gilbert3dpp-c

