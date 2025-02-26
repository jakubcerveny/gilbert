
CC := gcc
RUSTC := rustc
CFLAGS :=
OPT := -O3
RUSTOPT := -O

all: gilbert-c gilbert-rs

gilbert-c: gilbert.c
	$(CC) gilbert.c -o gilbert-c $(CFLAGS) $(OPT)

gilbert-rs: gilbert.rs
	$(RUSTC) --cfg 'feature="gilbert_main"' -o gilbert-rs gilbert.rs $(RUSTOPT)

.PHONY: clean all

clean:
	rm -f gilbert-c gilbert-rs
