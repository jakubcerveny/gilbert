#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-2-Clause
# Copyright (c) 2018 Jakub Červený


def gilbert3d(width, height, depth):
    """
    Generalized Hilbert ('Gilbert') space-filling curve for arbitrary-sized
    3D rectangular grids. Generates discrete 3D coordinates to fill a cuboid
    of size (width x height x depth). Even sizes are recommended in 3D.
    """

    if width >= height and width >= depth:
       yield from generate3d(0, 0, 0,
                             width, 0, 0,
                             0, height, 0,
                             0, 0, depth)

    elif height >= width and height >= depth:
       yield from generate3d(0, 0, 0,
                             0, height, 0,
                             width, 0, 0,
                             0, 0, depth)

    else: # depth >= width and depth >= height
       yield from generate3d(0, 0, 0,
                             0, 0, depth,
                             width, 0, 0,
                             0, height, 0)


def sgn(x):
    return -1 if x < 0 else (1 if x > 0 else 0)


def generate3d(x, y, z,
               ax, ay, az,
               bx, by, bz,
               cx, cy, cz):

    w = abs(ax + ay + az)
    h = abs(bx + by + bz)
    d = abs(cx + cy + cz)

    (dax, day, daz) = (sgn(ax), sgn(ay), sgn(az)) # unit major direction ("right")
    (dbx, dby, dbz) = (sgn(bx), sgn(by), sgn(bz)) # unit ortho direction ("forward")
    (dcx, dcy, dcz) = (sgn(cx), sgn(cy), sgn(cz)) # unit ortho direction ("up")

    # trivial row/column fills
    if h == 1 and d == 1:
        for i in range(0, w):
            yield(x, y, z)
            (x, y, z) = (x + dax, y + day, z + daz)
        return

    if w == 1 and d == 1:
        for i in range(0, h):
            yield(x, y, z)
            (x, y, z) = (x + dbx, y + dby, z + dbz)
        return

    if w == 1 and h == 1:
        for i in range(0, d):
            yield(x, y, z)
            (x, y, z) = (x + dcx, y + dcy, z + dcz)
        return

    (ax2, ay2, az2) = (ax//2, ay//2, az//2)
    (bx2, by2, bz2) = (bx//2, by//2, bz//2)
    (cx2, cy2, cz2) = (cx//2, cy//2, cz//2)

    w2 = abs(ax2 + ay2 + az2)
    h2 = abs(bx2 + by2 + bz2)
    d2 = abs(cx2 + cy2 + cz2)

    # prefer even steps
    if (w2 % 2) and (w > 2):
       (ax2, ay2, az2) = (ax2 + dax, ay2 + day, az2 + daz)

    if (h2 % 2) and (h > 2):
       (bx2, by2, bz2) = (bx2 + dbx, by2 + dby, bz2 + dbz)

    if (d2 % 2) and (d > 2):
       (cx2, cy2, cz2) = (cx2 + dcx, cy2 + dcy, cz2 + dcz)

    # wide case, split in w only
    if (2*w > 3*h) and (2*w > 3*d):
       yield from generate3d(x, y, z,
                             ax2, ay2, az2,
                             bx, by, bz,
                             cx, cy, cz)

       yield from generate3d(x+ax2, y+ay2, z+az2,
                             ax-ax2, ay-ay2, az-az2,
                             bx, by, bz,
                             cx, cy, cz)

    # do not split in d
    elif 3*h > 4*d:
       yield from generate3d(x, y, z,
                             bx2, by2, bz2,
                             cx, cy, cz,
                             ax2, ay2, az2)

       yield from generate3d(x+bx2, y+by2, z+bz2,
                             ax, ay, az,
                             bx-bx2, by-by2, bz-bz2,
                             cx, cy, cz)

       yield from generate3d(x+(ax-dax)+(bx2-dbx),
                             y+(ay-day)+(by2-dby),
                             z+(az-daz)+(bz2-dbz),
                             -bx2, -by2, -bz2,
                             cx, cy, cz,
                             -(ax-ax2), -(ay-ay2), -(az-az2))

    # do not split in h
    elif 3*d > 4*h:
       yield from generate3d(x, y, z,
                             cx2, cy2, cz2,
                             ax2, ay2, az2,
                             bx, by, bz)

       yield from generate3d(x+cx2, y+cy2, z+cz2,
                             ax, ay, az,
                             bx, by, bz,
                             cx-cx2, cy-cy2, cz-cz2)

       yield from generate3d(x+(ax-dax)+(cx2-dcx),
                             y+(ay-day)+(cy2-dcy),
                             z+(az-daz)+(cz2-dcz),
                             -cx2, -cy2, -cz2,
                             -(ax-ax2), -(ay-ay2), -(az-az2),
                             bx, by, bz)

    # regular case, split in all w/h/d
    else:
       yield from generate3d(x, y, z,
                             bx2, by2, bz2,
                             cx2, cy2, cz2,
                             ax2, ay2, az2)

       yield from generate3d(x+bx2, y+by2, z+bz2,
                             cx, cy, cz,
                             ax2, ay2, az2,
                             bx-bx2, by-by2, bz-bz2)

       yield from generate3d(x+(bx2-dbx)+(cx-dcx),
                             y+(by2-dby)+(cy-dcy),
                             z+(bz2-dbz)+(cz-dcz),
                             ax, ay, az,
                             -bx2, -by2, -bz2,
                             -(cx-cx2), -(cy-cy2), -(cz-cz2))

       yield from generate3d(x+(ax-dax)+bx2+(cx-dcx),
                             y+(ay-day)+by2+(cy-dcy),
                             z+(az-daz)+bz2+(cz-dcz),
                             -cx, -cy, -cz,
                             -(ax-ax2), -(ay-ay2), -(az-az2),
                             bx-bx2, by-by2, bz-bz2)

       yield from generate3d(x+(ax-dax)+(bx2-dbx),
                             y+(ay-day)+(by2-dby),
                             z+(az-daz)+(bz2-dbz),
                             -bx2, -by2, -bz2,
                             cx2, cy2, cz2,
                             -(ax-ax2), -(ay-ay2), -(az-az2))


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('width', type=int)
    parser.add_argument('height', type=int)
    parser.add_argument('depth', type=int)
    args = parser.parse_args()

    for x, y, z in gilbert3d(args.width, args.height, args.depth):
        print(x, y, z)
