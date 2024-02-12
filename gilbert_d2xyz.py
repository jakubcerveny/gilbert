#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-2-Clause
# Copyright (c) 2024 abetusk

def gilbert_d2xyz(idx, width, height, depth):
    """
    Generalized Hilbert ('Gilbert') space-filling curve for arbitrary-sized
    3D rectangular grids. Generates discrete 3D coordinates to fill a cuboid
    of size (width x height x depth). Even sizes are recommended in 3D.
    """

    if width >= height and width >= depth:
       return gilbert_d2xyz_r(idx, 0,
                              0, 0, 0,
                              width, 0, 0,
                              0, height, 0,
                              0, 0, depth)

    elif height >= width and height >= depth:
       return gilbert_d2xyz_r(idx, 0,
                              0, 0, 0,
                              0, height, 0,
                              width, 0, 0,
                              0, 0, depth)

    else: # depth >= width and depth >= height
       return gilbert_d2xyz_r(idx, 0,
                              0, 0, 0,
                              0, 0, depth,
                              width, 0, 0,
                              0, height, 0)

def sgn(x):
    return -1 if x < 0 else (1 if x > 0 else 0)

def gilbert_d2xyz_r(dst_idx, cur_idx,
                    x, y, z,
                    ax, ay, az,
                    bx, by, bz,
                    cx, cy, cz):

    w = abs(ax + ay + az)
    h = abs(bx + by + bz)
    d = abs(cx + cy + cz)

    (dax, day, daz) = (sgn(ax), sgn(ay), sgn(az)) # unit major direction ("right")
    (dbx, dby, dbz) = (sgn(bx), sgn(by), sgn(bz)) # unit ortho direction ("forward")
    (dcx, dcy, dcz) = (sgn(cx), sgn(cy), sgn(cz)) # unit ortho direction ("up")

    _dx = dax + dbx + dcx
    _dy = day + dby + dcy
    _dz = daz + dbz + dcz
    _di = dst_idx - cur_idx

    # trivial row/column fills
    if h == 1 and d == 1:
        return (x + dax*_di, y + day*_di, z + daz*_di)

    if w == 1 and d == 1:
        return (x + dbx*_di, y + dby*_di, z + dbz*_di)

    if w == 1 and h == 1:
        return (x + dcx*_di, y + dcy*_di, z + dcz*_di)

    (ax2, ay2, az2) = (ax//2, ay//2, az//2)
    (bx2, by2, bz2) = (bx//2, by//2, bz//2)
    (cx2, cy2, cz2) = (cx//2, cy//2, cz//2)

    w2 = abs(ax2 + ay2 + az2)
    h2 = abs(bx2 + by2 + bz2)
    d2 = abs(cx2 + cy2 + cz2)

    # prefer even steps
    if (w2 % 2) and (w > 2): (ax2, ay2, az2) = (ax2 + dax, ay2 + day, az2 + daz)
    if (h2 % 2) and (h > 2): (bx2, by2, bz2) = (bx2 + dbx, by2 + dby, bz2 + dbz)
    if (d2 % 2) and (d > 2): (cx2, cy2, cz2) = (cx2 + dcx, cy2 + dcy, cz2 + dcz)

    # wide case, split in w only
    if (2*w > 3*h) and (2*w > 3*d):
        nxt_idx = cur_idx + abs( (ax2 + ay2 + az2)*(bx + by + bz)*(cx + cy + cz) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x, y, z,
                                   ax2, ay2, az2,
                                   bx, by, bz,
                                   cx, cy, cz)
        cur_idx = nxt_idx

        return gilbert_d2xyz_r(dst_idx,cur_idx,
                               x+ax2, y+ay2, z+az2,
                               ax-ax2, ay-ay2, az-az2,
                               bx, by, bz,
                               cx, cy, cz)

    # do not split in d
    elif 3*h > 4*d:
        nxt_idx = cur_idx + abs( (bx2 + by2 + bz2)*(cx + cy + cz)*(ax2 + ay2 + az2) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x, y, z,
                                   bx2, by2, bz2,
                                   cx, cy, cz,
                                   ax2, ay2, az2)
        cur_idx = nxt_idx

        nxt_idx = cur_idx + abs( (ax + ay + az)*((bx - bx2) + (by - by2) + (bz - bz2))*(cx + cy + cz) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x+bx2, y+by2, z+bz2,
                                   ax, ay, az,
                                   bx-bx2, by-by2, bz-bz2,
                                   cx, cy, cz)
        cur_idx = nxt_idx

        return gilbert_d2xyz_r(dst_idx,cur_idx,
                               x+(ax-dax)+(bx2-dbx),
                               y+(ay-day)+(by2-dby),
                               z+(az-daz)+(bz2-dbz),
                               -bx2, -by2, -bz2,
                               cx, cy, cz,
                               -(ax-ax2), -(ay-ay2), -(az-az2))

    # do not split in h
    elif 3*d > 4*h:
        nxt_idx = cur_idx + abs( (cx2 + cy2 + cz2)*(ax2 + ay2 + az2)*(bx + by + bz) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x, y, z,
                                   cx2, cy2, cz2,
                                   ax2, ay2, az2,
                                   bx, by, bz)
        cur_idx = nxt_idx

        nxt_idx = cur_idx + abs( (ax + ay + az)*(bx + by + bz)*((cx - cx2) + (cy - cy2) + (cz - cz2)) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x+cx2, y+cy2, z+cz2,
                                   ax, ay, az,
                                   bx, by, bz,
                                   cx-cx2, cy-cy2, cz-cz2)
        cur_idx = nxt_idx

        return gilbert_d2xyz_r(dst_idx,cur_idx,
                               x+(ax-dax)+(cx2-dcx),
                               y+(ay-day)+(cy2-dcy),
                               z+(az-daz)+(cz2-dcz),
                               -cx2, -cy2, -cz2,
                               -(ax-ax2), -(ay-ay2), -(az-az2),
                               bx, by, bz)

    # regular case, split in all w/h/d
    else:
        nxt_idx = cur_idx + abs( (bx2 + by2 + bz2)*(cx2 + cy2 + cz2)*(ax2 + ay2 + az2) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x, y, z,
                                   bx2, by2, bz2,
                                   cx2, cy2, cz2,
                                   ax2, ay2, az2)
        cur_idx = nxt_idx

        nxt_idx = cur_idx + abs( (cx + cy + cz)*(ax2 + ay2 + az2)*((bx - bx2) + (by - by2) + (bz - bz2)) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x+bx2, y+by2, z+bz2,
                                   cx, cy, cz,
                                   ax2, ay2, az2,
                                   bx-bx2, by-by2, bz-bz2)
        cur_idx = nxt_idx

        nxt_idx = cur_idx + abs( (ax + ay + az)*(-bx2 - by2 - bz2)*(-(cx - cx2) - (cy - cy2) - (cz - cz2)) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx, cur_idx,
                                   x+(bx2-dbx)+(cx-dcx),
                                   y+(by2-dby)+(cy-dcy),
                                   z+(bz2-dbz)+(cz-dcz),
                                   ax, ay, az,
                                   -bx2, -by2, -bz2,
                                   -(cx-cx2), -(cy-cy2), -(cz-cz2))
        cur_idx = nxt_idx

        nxt_idx = cur_idx + abs( (-cx - cy - cz)*(-(ax - ax2) - (ay - ay2) - (az - az2))*((bx - bx2) + (by - by2) + (bz - bz2)) )
        if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
            return gilbert_d2xyz_r(dst_idx,cur_idx,
                                   x+(ax-dax)+bx2+(cx-dcx),
                                   y+(ay-day)+by2+(cy-dcy),
                                   z+(az-daz)+bz2+(cz-dcz),
                                   -cx, -cy, -cz,
                                   -(ax-ax2), -(ay-ay2), -(az-az2),
                                   bx-bx2, by-by2, bz-bz2)
        cur_idx = nxt_idx

        return gilbert_d2xyz_r(dst_idx,cur_idx,
                               x+(ax-dax)+(bx2-dbx),
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

    w = args.width
    h = args.height
    d = args.depth

    n = w*h*d

    for idx in range(n):
        (x,y,z) = gilbert_d2xyz(idx,w,h,d)
        print(x,y,z)

