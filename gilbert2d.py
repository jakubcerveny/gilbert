#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-2-Clause
# Copyright (c) 2018 Jakub Červený

import sys

def gilbert2d(width, height):
    """
    Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
    2D rectangular grids. Generates discrete 2D coordinates to fill a rectangle
    of size (width x height).
    """

    if width >= height:
        yield from generate2d(0, 0, width, 0, 0, height)
    else:
        yield from generate2d(0, 0, 0, height, width, 0)

def gilbertxy2d(x,y,w,h):
    """
    Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
    2D rectangular grids. Takes a discrete 2D coordinate and maps it to the
    index position on the gilbert curve.
    """

    return gilbertxy2d_r(0, x,y, 0,0, w,0, 0,h)

def gilbertd2xy(idx,w,h):
    """
    Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
    2D rectangular grids. Takes a position along the gilbert curve and returns
    its 2D (x,y) coordinate.
    """

    return gilbertd2xy_r(idx,0, 0,0, w,0, 0,h)


def sgn(x):
    return -1 if x < 0 else (1 if x > 0 else 0)

def inbounds(x,y, x_s,y_s, ax,ay, bx,by):

    dx = ax+bx
    dy = ay+by

    if (dx<0:
        if (x>x_s) or (x<=(x_s+dx)): return False
    else:
        if (x<x_s) or (x>=(x_s+dx)): return False

    if dy<0:
        if (y>y_s) or (y<=(y_s+dy)): return False
    else:
        if (y<y_s) or (y>=(y_s+dy)): return False

    return True


def gilbertd2xy_r( dst_idx, cur_idx, x,y, ax,ay, bx,by ):

    w = abs(ax + ay)
    h = abs(bx + by)

    (dax, day) = (sgn(ax), sgn(ay)) # unit major direction
    (dbx, dby) = (sgn(bx), sgn(by)) # unit orthogonal direction

    dx = dax+dbx
    dy = day+dby
    di = dst_idx - cur_idx

    if h == 1: return (x + dax*di, y + day*di)
    if w == 1: return (x + dbx*di, y + dby*di)

    (ax2, ay2) = (ax//2, ay//2)
    (bx2, by2) = (bx//2, by//2)

    w2 = abs(ax2 + ay2)
    h2 = abs(bx2 + by2)

    if 2*w > 3*h:
        if (w2 % 2) and (w > 2):
            # prefer even steps
            (ax2, ay2) = (ax2 + dax, ay2 + day)

        # long case: split in two parts only
        nxt_idx = cur_idx + abs((ax2 + ay2)*(bx + by))
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbertd2xy_r(dst_idx, cur_idx,  x, y, ax2, ay2, bx, by)
        cur_idx = nxt_idx

        return gilbertd2xy_r(dst_idx, cur_idx, x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by)

    else:
        if (h2 % 2) and (h > 2):
            # prefer even steps
            (bx2, by2) = (bx2 + dbx, by2 + dby)

        # standard case: one step up, one long horizontal, one step down
        nxt_idx = cur_idx + abs((bx2 + by2)*(ax2 + ay2))
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbertd2xy_r(dst_idx, cur_idx, x,y, bx2,by2, ax2,ay2)
        cur_idx = nxt_idx

        nxt_idx = cur_idx + abs((ax+ay)*((bx-bx2) + (by-by2)))
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbertd2xy_r(dst_idx, cur_idx, x+bx2,y+by2, ax,ay, bx-bx2,by-by2)
        cur_idx = nxt_idx

        return gilbertd2xy_r(dst_idx, cur_idx,
                       x+(ax-dax)+(bx2-dbx),
                       y+(ay-day)+(by2-dby),
                       -bx2, -by2,
                       -(ax-ax2), -(ay-ay2))

def gilbertxy2d_r( cur_idx, x_dst,y_dst, x,y, ax,ay, bx,by ):

    w = abs(ax + ay)
    h = abs(bx + by)

    (dax, day) = (sgn(ax), sgn(ay)) # unit major direction
    (dbx, dby) = (sgn(bx), sgn(by)) # unit orthogonal direction

    dx = dax+dbx
    dy = day+dby

    if h == 1:
        if (dax==0): return cur_idx + (dy*(y_dst-y))
        return cur_idx + (dx*(x_dst-x))

    if w == 1:
        if (dbx==0): return cur_idx + (dy*(y_dst-y))
        return cur_idx + (dx*(x_dst-x))

    (ax2, ay2) = (ax//2, ay//2)
    (bx2, by2) = (bx//2, by//2)

    w2 = abs(ax2 + ay2)
    h2 = abs(bx2 + by2)

    if 2*w > 3*h:
        if (w2 % 2) and (w > 2):
            # prefer even steps
            (ax2, ay2) = (ax2 + dax, ay2 + day)

        if inbounds( x_dst, y_dst, x,y, ax2,ay2, bx,by ):
            return gilbertxy2d_r(cur_idx, x_dst, y_dst, x, y, ax2, ay2, bx, by)

        cur_idx += abs((ax2 + ay2)*(bx + by))
        return gilbertxy2d_r(cur_idx, x_dst, y_dst, x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by)

    else:
        if (h2 % 2) and (h > 2):
            # prefer even steps
            (bx2, by2) = (bx2 + dbx, by2 + dby)

        # standard case: one step up, one long horizontal, one step down
        if inbounds( x_dst,y_dst, x,y, bx2,by2, ax2,ay2 ):
            return gilbertxy2d_r(cur_idx, x_dst,y_dst, x,y, bx2,by2, ax2,ay2)
        cur_idx += abs((bx2 + by2)*(ax2 + ay2))

        if inbounds( x_dst,y_dst, x+bx2, y+by2, ax, ay, bx-bx2, by-by2):
            return gilbertxy2d_r(cur_idx, x_dst,y_dst, x+bx2,y+by2, ax,ay, bx-bx2,by-by2)
        cur_idx += abs((ax+ay)*((bx-bx2) + (by-by2)))

        return gilbertxy2d_r(cur_idx, x_dst,y_dst,
                       x+(ax-dax)+(bx2-dbx),
                       y+(ay-day)+(by2-dby),
                       -bx2, -by2,
                       -(ax-ax2), -(ay-ay2))

def generate2d(x, y, ax, ay, bx, by):

    w = abs(ax + ay)
    h = abs(bx + by)

    (dax, day) = (sgn(ax), sgn(ay)) # unit major direction
    (dbx, dby) = (sgn(bx), sgn(by)) # unit orthogonal direction

    if h == 1:
        # trivial row fill
        for i in range(0, w):
            yield(x, y)
            (x, y) = (x + dax, y + day)
        return

    if w == 1:
        # trivial column fill
        for i in range(0, h):
            yield(x, y)
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
        yield from generate2d(x, y, ax2, ay2, bx, by)
        yield from generate2d(x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by)

    else:
        if (h2 % 2) and (h > 2):
            # prefer even steps
            (bx2, by2) = (bx2 + dbx, by2 + dby)

        # standard case: one step up, one long horizontal, one step down
        yield from generate2d(x, y, bx2, by2, ax2, ay2)
        yield from generate2d(x+bx2, y+by2, ax, ay, bx-bx2, by-by2)
        yield from generate2d(x+(ax-dax)+(bx2-dbx), y+(ay-day)+(by2-dby),
                              -bx2, -by2, -(ax-ax2), -(ay-ay2))


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('width', type=int)
    parser.add_argument('height', type=int)
    parser.add_argument('--op', type=str, required=False)
    args = parser.parse_args()

    width = args.width
    height = args.height
    op = args.op

    if op == "d2xy":

      for idx in range(width*height):
        (x,y) = gilbertd2xy(idx, width,height)
        print(x,y)

      sys.exit(0)

    elif op == "xy2d":

      for x in range(width):
        for y in range(height):
          idx = gilbertxy2d(x,y, width,height)
          print(idx,x,y)

      sys.exit(0)

    for x, y in gilbert2d(args.width, args.height):
        print(x, y)


