#!/usr/bin/env python3

import sys

def sgn(x):
    return -1 if x < 0 else (1 if x > 0 else 0)


def gilbert2d(x, y, ax, ay, bx, by):
    """
    Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
    2D rectangular grids.
    """

    w = abs(ax + ay)
    h = abs(bx + by)

    (dax, day) = (sgn(ax), sgn(ay)) # unit major direction
    (dbx, dby) = (sgn(bx), sgn(by)) # unit orthogonal direction

    if h == 1:
        # trivial row fill
        for i in range(0, w):
            print(x, y)
            (x, y) = (x + dax, y + day)
        return

    if w == 1:
        # trivial column fill
        for i in range(0, h):
            print(x, y)
            (x, y) = (x + dbx, y + dby)
        return

    (ax2, ay2) = (ax//2, ay//2)
    (bx2, by2) = (bx//2, by//2)

    w2 = abs(ax2 + ay2)
    h2 = abs(bx2 + by2)

    if 2*w > 3*h:
        if (w2 % 2) and (w > 2):
            # prefer even steps
            (ax2, ay2) = (ax2 + dax, ay2 + day)

        # long case: split in two parts only
        gilbert2d(x, y, ax2, ay2, bx, by)
        gilbert2d(x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by)

    else:
        if (h2 % 2) and (h > 2):
            # prefer even steps
            (bx2, by2) = (bx2 + dbx, by2 + dby)

        # standard case: one step up, one long horizontal, one step down
        gilbert2d(x, y, bx2, by2, ax2, ay2)
        gilbert2d(x+bx2, y+by2, ax, ay, bx-bx2, by-by2)
        gilbert2d(x+(ax-dax)+(bx2-dbx), y+(ay-day)+(by2-dby),
                 -bx2, -by2, -(ax-ax2), -(ay-ay2))


def main():
    width = int(sys.argv[1])
    height = int(sys.argv[2])

    if width >= height:
        gilbert2d(0, 0, width, 0, 0, height)
    else:
        gilbert2d(0, 0, 0, height, width, 0)


if __name__ == "__main__":
    main()
