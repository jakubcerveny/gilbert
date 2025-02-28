// SPDX-License-Identifier: BSD-2-Clause
// Copyright (c) 2025 наб <nabijaczleweli@nabijaczleweli.xyz>
//
// Suitable for inclusion, yielding mod gilbert.
// To build the test driver, use
// $ rustc --cfg 'feature="gilbert_main"' gilbert.rs

pub mod gilbert {
    pub fn xy2d(xy: (i32, i32), (w, h): (i32, i32)) -> i32 {
        if w >= h {
            xy2d_impl(0, xy, (0, 0), (w, 0), (0, h))
        } else {
            xy2d_impl(0, xy, (0, 0), (0, h), (w, 0))
        }
    }

    pub fn d2xy(idx: i32, (w, h): (i32, i32)) -> (i32, i32) {
        if w >= h {
            d2xy_impl(idx, 0, (0, 0), (w, 0), (0, h))
        } else {
            d2xy_impl(idx, 0, (0, 0), (0, h), (w, 0))
        }
    }

    pub fn xyz2d(xyz: (i32, i32, i32), (width, height, depth): (i32, i32, i32)) -> i32 {
        if (width >= height) && (width >= depth) {
            xyz2d_impl(0, xyz, (0, 0, 0), (width, 0, 0), (0, height, 0), (0, 0, depth))
        } else if (height >= width) && (height >= depth) {
            xyz2d_impl(0, xyz, (0, 0, 0), (0, height, 0), (width, 0, 0), (0, 0, depth))
        } else {
            // depth >= width and depth >= height
            xyz2d_impl(0, xyz, (0, 0, 0), (0, 0, depth), (width, 0, 0), (0, height, 0))
        }
    }

    pub fn d2xyz(idx: i32, (width, height, depth): (i32, i32, i32)) -> (i32, i32, i32) {
        if (width >= height) && (width >= depth) {
            d2xyz_impl(idx, 0, (0, 0, 0), (width, 0, 0), (0, height, 0), (0, 0, depth))
        } else if (height >= width) && (height >= depth) {
            d2xyz_impl(idx, 0, (0, 0, 0), (0, height, 0), (width, 0, 0), (0, 0, depth))
        } else {
            // depth >= width and depth >= height
            d2xyz_impl(idx, 0, (0, 0, 0), (0, 0, depth), (width, 0, 0), (0, height, 0))
        }
    }


    fn in_bounds2((x, y): (i32, i32), (x_s, y_s): (i32, i32), (ax, ay): (i32, i32), (bx, by): (i32, i32)) -> bool {
        let dx = ax + bx;
        let dy = ay + by;

        if dx < 0 {
            if (x > x_s) || (x <= (x_s + dx)) {
                return false;
            }
        } else {
            if (x < x_s) || (x >= (x_s + dx)) {
                return false;
            }
        }

        if dy < 0 {
            if (y > y_s) || (y <= (y_s + dy)) {
                return false;
            }
        } else {
            if (y < y_s) || (y >= (y_s + dy)) {
                return false;
            }
        }

        return true;
    }

    fn in_bounds3((x, y, z): (i32, i32, i32), (x_s, y_s, z_s): (i32, i32, i32), (ax, ay, az): (i32, i32, i32), (bx, by, bz): (i32, i32, i32),
                  (cx, cy, cz): (i32, i32, i32))
                  -> bool {
        let dx = ax + bx + cx;
        let dy = ay + by + cy;
        let dz = az + bz + cz;

        if dx < 0 {
            if (x > x_s) || (x <= (x_s + dx)) {
                return false;
            }
        } else {
            if (x < x_s) || (x >= (x_s + dx)) {
                return false;
            }
        }

        if dy < 0 {
            if (y > y_s) || (y <= (y_s + dy)) {
                return false;
            }
        } else {
            if (y < y_s) || (y >= (y_s + dy)) {
                return false;
            }
        }

        if dz < 0 {
            if (z > z_s) || (z <= (z_s + dz)) {
                return false;
            }
        } else {
            if (z < z_s) || (z >= (z_s + dz)) {
                return false;
            }
        }

        return true;
    }



    fn d2xy_impl(dst_idx: i32, mut cur_idx: i32, (x, y): (i32, i32), (ax, ay): (i32, i32), (bx, by): (i32, i32)) -> (i32, i32) {
        let w = (ax + ay).abs();
        let h = (bx + by).abs();

        let (dax, day) = (ax.signum(), ay.signum()); // unit major direction
        let (dbx, dby) = (bx.signum(), by.signum()); // unit orthogonal direction

        let di = dst_idx - cur_idx;

        if h == 1 {
            return (x + dax * di, y + day * di);
        }

        if w == 1 {
            return (x + dbx * di, y + dby * di);
        }

        // floor function
        let (mut ax2, mut ay2) = (ax >> 1, ay >> 1);
        let (mut bx2, mut by2) = (bx >> 1, by >> 1);

        let w2 = (ax2 + ay2).abs();
        let h2 = (bx2 + by2).abs();

        if (2 * w) > (3 * h) {
            if (w2 & 1) != 0 && (w > 2) {
                // prefer even steps
                ax2 += dax;
                ay2 += day;
            }

            // long case: split in two parts only
            let nxt_idx = cur_idx + ((ax2 + ay2) * (bx + by)).abs();
            if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
                return d2xy_impl(dst_idx, cur_idx, (x, y), (ax2, ay2), (bx, by));
            }
            cur_idx = nxt_idx;

            return d2xy_impl(dst_idx, cur_idx, (x + ax2, y + ay2), (ax - ax2, ay - ay2), (bx, by));
        }

        if (h2 & 1) != 0 && (h > 2) {
            // prefer even steps
            bx2 += dbx;
            by2 += dby;
        }

        // standard case: one step up, one long horizontal, one step down
        let nxt_idx = cur_idx + ((bx2 + by2) * (ax2 + ay2)).abs();
        if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
            return d2xy_impl(dst_idx, cur_idx, (x, y), (bx2, by2), (ax2, ay2));
        }
        cur_idx = nxt_idx;

        let nxt_idx = cur_idx + ((ax + ay) * ((bx - bx2) + (by - by2))).abs();
        if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
            return d2xy_impl(dst_idx, cur_idx, (x + bx2, y + by2), (ax, ay), (bx - bx2, by - by2));
        }

        cur_idx = nxt_idx;

        return d2xy_impl(dst_idx,
                         cur_idx,
                         (x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby)),
                         (-bx2, -by2),
                         (-(ax - ax2), -(ay - ay2)));
    }

    fn xy2d_impl(mut cur_idx: i32, (x_dst, y_dst): (i32, i32), (x, y): (i32, i32), (ax, ay): (i32, i32), (bx, by): (i32, i32)) -> i32 {
        let w = (ax + ay).abs();
        let h = (bx + by).abs();

        let (dax, day) = (ax.signum(), ay.signum()); // unit major direction
        let (dbx, dby) = (bx.signum(), by.signum()); // unit orthogonal direction

        let dx = dax + dbx;
        let dy = day + dby;

        if h == 1 {
            if dax == 0 {
                return cur_idx + (dy * (y_dst - y));
            } else {
                return cur_idx + (dx * (x_dst - x));
            }
        }

        if w == 1 {
            if dbx == 0 {
                return cur_idx + (dy * (y_dst - y));
            } else {
                return cur_idx + (dx * (x_dst - x));
            }
        }

        let (mut ax2, mut ay2) = (ax >> 1, ay >> 1);
        let (mut bx2, mut by2) = (bx >> 1, by >> 1);

        let w2 = (ax2 + ay2).abs();
        let h2 = (bx2 + by2).abs();

        if (2 * w) > (3 * h) {
            if (w2 & 1) != 0 && (w > 2) {
                // prefer even steps
                ax2 += dax;
                ay2 += day;
            }

            if in_bounds2((x_dst, y_dst), (x, y), (ax2, ay2), (bx, by)) {
                return xy2d_impl(cur_idx, (x_dst, y_dst), (x, y), (ax2, ay2), (bx, by));
            }
            cur_idx += ((ax2 + ay2) * (bx + by)).abs();

            return xy2d_impl(cur_idx, (x_dst, y_dst), (x + ax2, y + ay2), (ax - ax2, ay - ay2), (bx, by));
        }

        if (h2 & 1) != 0 && (h > 2) {
            // prefer even steps
            bx2 += dbx;
            by2 += dby;
        }

        // standard case: one step up, one long horizontal, one step down
        if in_bounds2((x_dst, y_dst), (x, y), (bx2, by2), (ax2, ay2)) {
            return xy2d_impl(cur_idx, (x_dst, y_dst), (x, y), (bx2, by2), (ax2, ay2));
        }
        cur_idx += ((bx2 + by2) * (ax2 + ay2)).abs();

        if in_bounds2((x_dst, y_dst), (x + bx2, y + by2), (ax, ay), (bx - bx2, by - by2)) {
            return xy2d_impl(cur_idx, (x_dst, y_dst), (x + bx2, y + by2), (ax, ay), (bx - bx2, by - by2));
        }
        cur_idx += ((ax + ay) * ((bx - bx2) + (by - by2))).abs();

        return xy2d_impl(cur_idx,
                         (x_dst, y_dst),
                         (x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby)),
                         (-bx2, -by2),
                         (-(ax - ax2), -(ay - ay2)));
    }



    fn d2xyz_impl(dst_idx: i32, mut cur_idx: i32, (x, y, z): (i32, i32, i32), (ax, ay, az): (i32, i32, i32), (bx, by, bz): (i32, i32, i32),
                  (cx, cy, cz): (i32, i32, i32))
                  -> (i32, i32, i32) {
        let w = (ax + ay + az).abs();
        let h = (bx + by + bz).abs();
        let d = (cx + cy + cz).abs();

        let (dax, day, daz) = (ax.signum(), ay.signum(), az.signum()); // unit major direction "right"
        let (dbx, dby, dbz) = (bx.signum(), by.signum(), bz.signum()); // unit ortho direction "forward"
        let (dcx, dcy, dcz) = (cx.signum(), cy.signum(), cz.signum()); // unit ortho direction "up"

        let _di = dst_idx - cur_idx;

        // trivial row/column fills
        if (h == 1) && (d == 1) {
            return (x + dax * _di, y + day * _di, z + daz * _di);
        }

        if (w == 1) && (d == 1) {
            return (x + dbx * _di, y + dby * _di, z + dbz * _di);
        }

        if (w == 1) && (h == 1) {
            return (x + dcx * _di, y + dcy * _di, z + dcz * _di);
        }

        let (mut ax2, mut ay2, mut az2) = (ax >> 1, ay >> 1, az >> 1);
        let (mut bx2, mut by2, mut bz2) = (bx >> 1, by >> 1, bz >> 1);
        let (mut cx2, mut cy2, mut cz2) = (cx >> 1, cy >> 1, cz >> 1);

        let w2 = (ax2 + ay2 + az2).abs();
        let h2 = (bx2 + by2 + bz2).abs();
        let d2 = (cx2 + cy2 + cz2).abs();

        // prefer even steps
        if (w2 & 1) != 0 && (w > 2) {
            ax2 += dax;
            ay2 += day;
            az2 += daz;
        }
        if (h2 & 1) != 0 && (h > 2) {
            bx2 += dbx;
            by2 += dby;
            bz2 += dbz;
        }
        if (d2 & 1) != 0 && (d > 2) {
            cx2 += dcx;
            cy2 += dcy;
            cz2 += dcz;
        }

        // wide case, split in w only
        if ((2 * w) > (3 * h)) && ((2 * w) > (3 * d)) {
            let nxt_idx = cur_idx + ((ax2 + ay2 + az2) * (bx + by + bz) * (cx + cy + cz)).abs();
            if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
                return d2xyz_impl(dst_idx, cur_idx, (x, y, z), (ax2, ay2, az2), (bx, by, bz), (cx, cy, cz));
            }
            cur_idx = nxt_idx;


            return d2xyz_impl(dst_idx,
                              cur_idx,
                              (x + ax2, y + ay2, z + az2),
                              (ax - ax2, ay - ay2, az - az2),
                              (bx, by, bz),
                              (cx, cy, cz));
        }
        // do not split in d
        else if (3 * h) > (4 * d) {
            let nxt_idx = cur_idx + ((bx2 + by2 + bz2) * (cx + cy + cz) * (ax2 + ay2 + az2)).abs();
            if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
                return d2xyz_impl(dst_idx, cur_idx, (x, y, z), (bx2, by2, bz2), (cx, cy, cz), (ax2, ay2, az2));
            }
            cur_idx = nxt_idx;

            let nxt_idx = cur_idx + ((ax + ay + az) * ((bx - bx2) + (by - by2) + (bz - bz2)) * (cx + cy + cz)).abs();
            if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
                return d2xyz_impl(dst_idx,
                                  cur_idx,
                                  (x + bx2, y + by2, z + bz2),
                                  (ax, ay, az),
                                  (bx - bx2, by - by2, bz - bz2),
                                  (cx, cy, cz));
            }
            cur_idx = nxt_idx;

            return d2xyz_impl(dst_idx,
                              cur_idx,
                              (x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), z + (az - daz) + (bz2 - dbz)),
                              (-bx2, -by2, -bz2),
                              (cx, cy, cz),
                              (-(ax - ax2), -(ay - ay2), -(az - az2)));
        }
        // do not split in h
        else if (3 * d) > (4 * h) {
            let nxt_idx = cur_idx + ((cx2 + cy2 + cz2) * (ax2 + ay2 + az2) * (bx + by + bz)).abs();
            if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
                return d2xyz_impl(dst_idx, cur_idx, (x, y, z), (cx2, cy2, cz2), (ax2, ay2, az2), (bx, by, bz));
            }
            cur_idx = nxt_idx;

            let nxt_idx = cur_idx + ((ax + ay + az) * (bx + by + bz) * ((cx - cx2) + (cy - cy2) + (cz - cz2))).abs();
            if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
                return d2xyz_impl(dst_idx,
                                  cur_idx,
                                  (x + cx2, y + cy2, z + cz2),
                                  (ax, ay, az),
                                  (bx, by, bz),
                                  (cx - cx2, cy - cy2, cz - cz2));
            }
            cur_idx = nxt_idx;

            return d2xyz_impl(dst_idx,
                              cur_idx,
                              (x + (ax - dax) + (cx2 - dcx), y + (ay - day) + (cy2 - dcy), z + (az - daz) + (cz2 - dcz)),
                              (-cx2, -cy2, -cz2),
                              (-(ax - ax2), -(ay - ay2), -(az - az2)),
                              (bx, by, bz));

        }

        // regular case, split in all w/h/d
        let nxt_idx = cur_idx + ((bx2 + by2 + bz2) * (cx2 + cy2 + cz2) * (ax2 + ay2 + az2)).abs();
        if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
            return d2xyz_impl(dst_idx, cur_idx, (x, y, z), (bx2, by2, bz2), (cx2, cy2, cz2), (ax2, ay2, az2));
        }
        cur_idx = nxt_idx;

        let nxt_idx = cur_idx + ((cx + cy + cz) * (ax2 + ay2 + az2) * ((bx - bx2) + (by - by2) + (bz - bz2))).abs();
        if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
            return d2xyz_impl(dst_idx,
                              cur_idx,
                              (x + bx2, y + by2, z + bz2),
                              (cx, cy, cz),
                              (ax2, ay2, az2),
                              (bx - bx2, by - by2, bz - bz2));
        }
        cur_idx = nxt_idx;

        let nxt_idx = cur_idx + ((ax + ay + az) * (-bx2 - by2 - bz2) * (-(cx - cx2) - (cy - cy2) - (cz - cz2))).abs();
        if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
            return d2xyz_impl(dst_idx,
                              cur_idx,
                              (x + (bx2 - dbx) + (cx - dcx), y + (by2 - dby) + (cy - dcy), z + (bz2 - dbz) + (cz - dcz)),
                              (ax, ay, az),
                              (-bx2, -by2, -bz2),
                              (-(cx - cx2), -(cy - cy2), -(cz - cz2)));
        }
        cur_idx = nxt_idx;

        let nxt_idx = cur_idx + ((-cx - cy - cz) * (-(ax - ax2) - (ay - ay2) - (az - az2)) * ((bx - bx2) + (by - by2) + (bz - bz2))).abs();
        if (cur_idx <= dst_idx) && (dst_idx < nxt_idx) {
            return d2xyz_impl(dst_idx,
                              cur_idx,
                              (x + (ax - dax) + bx2 + (cx - dcx), y + (ay - day) + by2 + (cy - dcy), z + (az - daz) + bz2 + (cz - dcz)),
                              (-cx, -cy, -cz),
                              (-(ax - ax2), -(ay - ay2), -(az - az2)),
                              (bx - bx2, by - by2, bz - bz2));
        }
        cur_idx = nxt_idx;

        return d2xyz_impl(dst_idx,
                          cur_idx,
                          (x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), z + (az - daz) + (bz2 - dbz)),
                          (-bx2, -by2, -bz2),
                          (cx2, cy2, cz2),
                          (-(ax - ax2), -(ay - ay2), -(az - az2)));
    }



    fn xyz2d_impl(mut cur_idx: i32, (x_dst, y_dst, z_dst): (i32, i32, i32), (x, y, z): (i32, i32, i32), (ax, ay, az): (i32, i32, i32),
                  (bx, by, bz): (i32, i32, i32), (cx, cy, cz): (i32, i32, i32))
                  -> i32 {
        let w = (ax + ay + az).abs();
        let h = (bx + by + bz).abs();
        let d = (cx + cy + cz).abs();

        let (dax, day, daz) = (ax.signum(), ay.signum(), az.signum()); // unit major direction "right"
        let (dbx, dby, dbz) = (bx.signum(), by.signum(), bz.signum()); // unit ortho direction "forward"
        let (dcx, dcy, dcz) = (cx.signum(), cy.signum(), cz.signum()); // unit ortho direction "up"

        // trivial row/column fills
        if (h == 1) && (d == 1) {
            return cur_idx + (dax * (x_dst - x)) + (day * (y_dst - y)) + (daz * (z_dst - z));
        }

        if (w == 1) && (d == 1) {
            return cur_idx + (dbx * (x_dst - x)) + (dby * (y_dst - y)) + (dbz * (z_dst - z));
        }

        if (w == 1) && (h == 1) {
            return cur_idx + (dcx * (x_dst - x)) + (dcy * (y_dst - y)) + (dcz * (z_dst - z));
        }

        let (mut ax2, mut ay2, mut az2) = (ax >> 1, ay >> 1, az >> 1);
        let (mut bx2, mut by2, mut bz2) = (bx >> 1, by >> 1, bz >> 1);
        let (mut cx2, mut cy2, mut cz2) = (cx >> 1, cy >> 1, cz >> 1);

        let w2 = (ax2 + ay2 + az2).abs();
        let h2 = (bx2 + by2 + bz2).abs();
        let d2 = (cx2 + cy2 + cz2).abs();

        // prefer even steps
        if (w2 & 1) != 0 && (w > 2) {
            ax2 += dax;
            ay2 += day;
            az2 += daz;
        }

        if (h2 & 1) != 0 && (h > 2) {
            bx2 += dbx;
            by2 += dby;
            bz2 += dbz;
        }

        if (d2 & 1) != 0 && (d > 2) {
            cx2 += dcx;
            cy2 += dcy;
            cz2 += dcz;
        }

        // wide case, split in w only
        if (2 * w > 3 * h) && (2 * w > 3 * d) {
            if in_bounds3((x_dst, y_dst, z_dst), (x, y, z), (ax2, ay2, az2), (bx, by, bz), (cx, cy, cz)) {
                return xyz2d_impl(cur_idx, (x_dst, y_dst, z_dst), (x, y, z), (ax2, ay2, az2), (bx, by, bz), (cx, cy, cz));
            }
            cur_idx += ((ax2 + ay2 + az2) * (bx + by + bz) * (cx + cy + cz)).abs();

            return xyz2d_impl(cur_idx,
                              (x_dst, y_dst, z_dst),
                              (x + ax2, y + ay2, z + az2),
                              (ax - ax2, ay - ay2, az - az2),
                              (bx, by, bz),
                              (cx, cy, cz));
        }
        // do not split in d
        else if (3 * h) > (4 * d) {
            if in_bounds3((x_dst, y_dst, z_dst), (x, y, z), (bx2, by2, bz2), (cx, cy, cz), (ax2, ay2, az2)) {
                return xyz2d_impl(cur_idx, (x_dst, y_dst, z_dst), (x, y, z), (bx2, by2, bz2), (cx, cy, cz), (ax2, ay2, az2));
            }
            cur_idx += ((bx2 + by2 + bz2) * (cx + cy + cz) * (ax2 + ay2 + az2)).abs();

            if in_bounds3((x_dst, y_dst, z_dst),
                          (x + bx2, y + by2, z + bz2),
                          (ax, ay, az),
                          (bx - bx2, by - by2, bz - bz2),
                          (cx, cy, cz)) {
                return xyz2d_impl(cur_idx,
                                  (x_dst, y_dst, z_dst),
                                  (x + bx2, y + by2, z + bz2),
                                  (ax, ay, az),
                                  (bx - bx2, by - by2, bz - bz2),
                                  (cx, cy, cz));
            }
            cur_idx += ((ax + ay + az) * ((bx - bx2) + (by - by2) + (bz - bz2)) * (cx + cy + cz)).abs();

            return xyz2d_impl(cur_idx,
                              (x_dst, y_dst, z_dst),
                              (x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), z + (az - daz) + (bz2 - dbz)),
                              (-bx2, -by2, -bz2),
                              (cx, cy, cz),
                              (-(ax - ax2), -(ay - ay2), -(az - az2)));
        }
        // do not split in h
        else if (3 * d) > (4 * h) {
            if in_bounds3((x_dst, y_dst, z_dst), (x, y, z), (cx2, cy2, cz2), (ax2, ay2, az2), (bx, by, bz)) {
                return xyz2d_impl(cur_idx, (x_dst, y_dst, z_dst), (x, y, z), (cx2, cy2, cz2), (ax2, ay2, az2), (bx, by, bz));
            }
            cur_idx += ((cx2 + cy2 + cz2) * (ax2 + ay2 + az2) * (bx + by + bz)).abs();

            if in_bounds3((x_dst, y_dst, z_dst),
                          (x + cx2, y + cy2, z + cz2),
                          (ax, ay, az),
                          (bx, by, bz),
                          (cx - cx2, cy - cy2, cz - cz2)) {
                return xyz2d_impl(cur_idx,
                                  (x_dst, y_dst, z_dst),
                                  (x + cx2, y + cy2, z + cz2),
                                  (ax, ay, az),
                                  (bx, by, bz),
                                  (cx - cx2, cy - cy2, cz - cz2));
            }
            cur_idx += ((ax + ay + az) * (bx + by + bz) * ((cx - cx2) + (cy - cy2) + (cz - cz2))).abs();

            return xyz2d_impl(cur_idx,
                              (x_dst, y_dst, z_dst),
                              (x + (ax - dax) + (cx2 - dcx), y + (ay - day) + (cy2 - dcy), z + (az - daz) + (cz2 - dcz)),
                              (-cx2, -cy2, -cz2),
                              (-(ax - ax2), -(ay - ay2), -(az - az2)),
                              (bx, by, bz));

        }

        // regular case, split in all w/h/d
        if in_bounds3((x_dst, y_dst, z_dst), (x, y, z), (bx2, by2, bz2), (cx2, cy2, cz2), (ax2, ay2, az2)) {
            return xyz2d_impl(cur_idx, (x_dst, y_dst, z_dst), (x, y, z), (bx2, by2, bz2), (cx2, cy2, cz2), (ax2, ay2, az2));
        }
        cur_idx += ((bx2 + by2 + bz2) * (cx2 + cy2 + cz2) * (ax2 + ay2 + az2)).abs();

        if in_bounds3((x_dst, y_dst, z_dst),
                      (x + bx2, y + by2, z + bz2),
                      (cx, cy, cz),
                      (ax2, ay2, az2),
                      (bx - bx2, by - by2, bz - bz2)) {
            return xyz2d_impl(cur_idx,
                              (x_dst, y_dst, z_dst),
                              (x + bx2, y + by2, z + bz2),
                              (cx, cy, cz),
                              (ax2, ay2, az2),
                              (bx - bx2, by - by2, bz - bz2));
        }
        cur_idx += ((cx + cy + cz) * (ax2 + ay2 + az2) * ((bx - bx2) + (by - by2) + (bz - bz2))).abs();

        if in_bounds3((x_dst, y_dst, z_dst),
                      (x + (bx2 - dbx) + (cx - dcx), y + (by2 - dby) + (cy - dcy), z + (bz2 - dbz) + (cz - dcz)),
                      (ax, ay, az),
                      (-bx2, -by2, -bz2),
                      (-(cx - cx2), -(cy - cy2), -(cz - cz2))) {
            return xyz2d_impl(cur_idx,
                              (x_dst, y_dst, z_dst),
                              (x + (bx2 - dbx) + (cx - dcx), y + (by2 - dby) + (cy - dcy), z + (bz2 - dbz) + (cz - dcz)),
                              (ax, ay, az),
                              (-bx2, -by2, -bz2),
                              (-(cx - cx2), -(cy - cy2), -(cz - cz2)));
        }
        cur_idx += ((ax + ay + az) * (-bx2 - by2 - bz2) * (-(cx - cx2) - (cy - cy2) - (cz - cz2))).abs();

        if in_bounds3((x_dst, y_dst, z_dst),
                      (x + (ax - dax) + bx2 + (cx - dcx), y + (ay - day) + by2 + (cy - dcy), z + (az - daz) + bz2 + (cz - dcz)),
                      (-cx, -cy, -cz),
                      (-(ax - ax2), -(ay - ay2), -(az - az2)),
                      (bx - bx2, by - by2, bz - bz2)) {
            return xyz2d_impl(cur_idx,
                              (x_dst, y_dst, z_dst),
                              (x + (ax - dax) + bx2 + (cx - dcx), y + (ay - day) + by2 + (cy - dcy), z + (az - daz) + bz2 + (cz - dcz)),
                              (-cx, -cy, -cz),
                              (-(ax - ax2), -(ay - ay2), -(az - az2)),
                              (bx - bx2, by - by2, bz - bz2));
        }
        cur_idx += ((-cx - cy - cz) * (-(ax - ax2) - (ay - ay2) - (az - az2)) * ((bx - bx2) + (by - by2) + (bz - bz2))).abs();

        return xyz2d_impl(cur_idx,
                          (x_dst, y_dst, z_dst),
                          (x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), z + (az - daz) + (bz2 - dbz)),
                          (-bx2, -by2, -bz2),
                          (cx2, cy2, cz2),
                          (-(ax - ax2), -(ay - ay2), -(az - az2)));
    }
}


#[cfg(feature="gilbert_main")]
fn main() {
    static USAGE: &str = "gilbert xy2d|d2xy|xyz2d|d2xyz width height [depth]\ndepth defaults to 1 in 3D\n";
    let mut args = std::env::args().skip(1);
    let op = args.next().expect(USAGE);
    let w = args.next().expect(USAGE).parse().unwrap();
    let h = args.next().expect(USAGE).parse().unwrap();
    let d = args.next().as_deref().unwrap_or("1").parse().unwrap();

    match &op[..] {
        "xy2d" => {
            for x in 0..w {
                for y in 0..h {
                    let idx = gilbert::xy2d((x, y), (w, h));
                    println!("{} {} {}", idx, x, y);
                }
            }
        }
        "d2xy" => {
            for idx in 0..(w * h) {
                let (x, y) = gilbert::d2xy(idx, (w, h));
                println!("{} {}", x, y);
            }
        }
        "xyz2d" => {
            for x in 0..w {
                for y in 0..h {
                    for z in 0..d {
                        let idx = gilbert::xyz2d((x, y, z), (w, h, d));
                        println!("{} {} {} {}", idx, x, y, z);
                    }
                }
            }
        }

        "d2xyz" => {
            for idx in 0..(w * h * d) {
                let (x, y, z) = gilbert::d2xyz(idx, (w, h, d));
                println!("{} {} {}", x, y, z);
            }
        }
        _ => panic!("{}", USAGE),
    }
}
