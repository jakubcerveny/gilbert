
CC := gcc
RUSTC := rustc
CFLAGS :=
OPT := -O3
RUSTOPT := -O

#all: gilbert-c gilbert-rs
all: gilbert3dpp-c

gilbert3dpp-c: gilbert3dpp.c
	$(CC) gilbert3dpp.c -o gilbert3dpp-c $(CFLAGS) $(OPT)

#gilbert-rs: gilbert.rs
#	$(RUSTC) --cfg 'feature="gilbert_main"' -o gilbert-rs gilbert.rs $(RUSTOPT)

.PHONY: clean all

clean:
	rm -f gilbert3dpp-c gilbert3dpp-rs

