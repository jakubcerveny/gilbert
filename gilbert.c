// SPDX-License-Identifier: BSD-2-Clause
// Copyright (c) 2018 Jakub Červený

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int gilbertd2xy_r(int dst_idx, int cur_idx,
                  int *xres, int *yres,
                  int ax,int ay,
                  int bx,int by );

int gilbertxy2d_r(int cur_idx,
                  int x_dst, int y_dst,
                  int x, int y,
                  int ax, int ay,
                  int bx,int by );

int gilbertxy2d(int x, int y, int w, int h) {
  if (w >= h) {
    return gilbertxy2d_r(0, x,y, 0,0, w,0, 0,h);
  }
  return gilbertxy2d_r(0, x,y, 0,0, 0,h, w,0);
}

int gilbertd2xy(int *x, int *y, int idx,int w,int h) {
  *x = 0;
  *y = 0;

  if (w >= h) {
    return gilbertd2xy_r(idx,0, x,y, w,0, 0,h);
  }
  return gilbertd2xy_r(idx,0, x,y, 0,h, w,0);
}

int gilbertd2xyz_r(int dst_idx, int cur_idx,
                   int *x, int *y, int *z,
                   int ax, int ay, int az,
                   int bx, int by, int bz,
                   int cx, int cy, int cz);

int gilbertxyz2d_r(int cur_idx,
                   int x_dst, int y_dst, int z_dst,
                   int x,  int y,  int z,
                   int ax, int ay, int az,
                   int bx, int by, int bz,
                   int cx, int cy, int cz);


int gilbertxyz2d(int x, int y, int z, int width, int height, int depth) {
  if ((width >= height) && (width >= depth)) {
    return gilbertxyz2d_r( 0,x,y,z,
                           0, 0, 0,
                           width, 0, 0,
                           0, height, 0,
                           0, 0, depth);
  }
  else if ((height >= width) && (height >= depth)) {
    return gilbertxyz2d_r( 0,x,y,z,
                           0, 0, 0,
                           0, height, 0,
                           width, 0, 0,
                           0, 0, depth);
  }

  // depth >= width and depth >= height
  return gilbertxyz2d_r( 0,x,y,z,
                         0, 0, 0,
                         0, 0, depth,
                         width, 0, 0,
                         0, height, 0);
}

int gilbertd2xyz(int *x, int *y, int *z, int idx, int width, int height, int depth) {

  *x = 0;
  *y = 0;
  *z = 0;

  if ((width >= height) && (width >= depth)) {
    return gilbertd2xyz_r( idx, 0,
                           x,y,z,
                           width, 0, 0,
                           0, height, 0,
                           0, 0, depth);
  }

  else if ((height >= width) && (height >= depth)) {
    return gilbertd2xyz_r( idx, 0,
                           x,y,z,
                           0, height, 0,
                           width, 0, 0,
                           0, 0, depth);
  }

  // depth >= width and depth >= height
  return gilbertd2xyz_r( idx, 0,
                         x,y,z,
                         0, 0, depth,
                         width, 0, 0,
                         0, height, 0);
}


int sgn(int x) {
  if (x<0) { return -1; }
  if (x>0) { return  1; }
  return 0;
}

int inbounds2(int x,  int y,
              int x_s,int y_s,
              int ax, int ay,
              int bx, int by) {
  int dx,dy;

  dx = ax+bx;
  dy = ay+by;

  if (dx<0) {
    if ((x>x_s) || (x<=(x_s+dx))) { return 0; }
  }
  else {
    if ((x<x_s) || (x>=(x_s+dx))) { return 0; }
  }

  if (dy<0) {
    if ((y>y_s) || (y<=(y_s+dy))) { return 0; }
  }
  else {
    if ((y<y_s) || (y>=(y_s+dy))) { return 0; }
  }

  return 1;
}

int inbounds3(int x,  int y,  int z,
              int x_s,int y_s,int z_s,
              int ax, int ay, int az,
              int bx, int by, int bz,
              int cx, int cy, int cz) {
  int dx, dy, dz;

  dx = ax+bx+cx;
  dy = ay+by+cy;
  dz = az+bz+cz;

  if (dx<0) {
    if ((x>x_s) || (x<=(x_s+dx))) { return 0; }
  }
  else {
    if ((x<x_s) || (x>=(x_s+dx))) { return 0; }
  }

  if (dy<0) {
    if ((y>y_s) || (y<=(y_s+dy))) { return 0; }
  }
  else {
    if ((y<y_s) || (y>=(y_s+dy))) { return 0; }
  }

  if (dz<0) {
    if ((z>z_s) || (z<=(z_s+dz))) { return 0; }
  }
  else {
    if ((z<z_s) || (z>=(z_s+dz))) { return 0; }
  }

  return 1;
}


int gilbertd2xy_r(int dst_idx, int cur_idx,
                  int *xres, int *yres,
                  int ax,int ay,
                  int bx,int by ) {
  int nxt_idx;
  int w,h, x,y,
      dax,day,
      dbx,dby,
      dx,dy,di;
  int ax2,ay2, bx2,by2, w2,h2;

  w = abs(ax + ay);
  h = abs(bx + by);

  x = *xres;
  y = *yres;

  // unit major direction
  dax = sgn(ax);
  day = sgn(ay);

  // unit orthogonal direction
  dbx = sgn(bx);
  dby = sgn(by);

  dx = dax+dbx;
  dy = day+dby;
  di = dst_idx - cur_idx;

  if (h == 1) {
    *xres = x + dax*di;
    *yres = y + day*di;
    return 0;
  }

  if (w == 1) {
    *xres = x + dbx*di;
    *yres = y + dby*di;
    return 0;
  }

  ax2 = ax/2; ay2 = ay/2;
  bx2 = bx/2; by2 = by/2;

  w2 = abs(ax2 + ay2);
  h2 = abs(bx2 + by2);

  if ((2*w) > (3*h)) {
    if ((w2 % 2) && (w > 2)) {
      // prefer even steps
      ax2 += dax;
      ay2 += day;
    }

    // long case: split in two parts only
    nxt_idx = cur_idx + abs((ax2 + ay2)*(bx + by));
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      *xres = x;
      *yres = y;
      return gilbertd2xy_r(dst_idx, cur_idx,  xres, yres, ax2, ay2, bx, by);
    }
    cur_idx = nxt_idx;

    *xres = x+ax2;
    *yres = y+ay2;
    return gilbertd2xy_r(dst_idx, cur_idx, xres, yres, ax-ax2, ay-ay2, bx, by);
  }

  if ((h2 % 2) && (h > 2)) {
    // prefer even steps
    bx2 += dbx;
    by2 += dby;
  }

  // standard case: one step up, one long horizontal, one step down
  nxt_idx = cur_idx + abs((bx2 + by2)*(ax2 + ay2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    *xres = x;
    *yres = y;
    return gilbertd2xy_r(dst_idx, cur_idx, xres,yres, bx2,by2, ax2,ay2);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + abs((ax+ay)*((bx-bx2) + (by-by2)));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    *xres = x+bx2;
    *yres = y+by2;
    return gilbertd2xy_r(dst_idx, cur_idx, xres,yres, ax,ay, bx-bx2,by-by2);
  }
  cur_idx = nxt_idx;

  *xres = x+(ax-dax)+(bx2-dbx);
  *yres = y+(ay-day)+(by2-dby);

  return gilbertd2xy_r( dst_idx, cur_idx,
                        xres,yres,
                        -bx2, -by2,
                        -(ax-ax2), -(ay-ay2));
}

int gilbertxy2d_r(int cur_idx,
                  int x_dst, int y_dst,
                  int x, int y,
                  int ax, int ay,
                  int bx,int by ) {
  int dax,day,  dbx,dby,
      ax2,ay2,  bx2,by2;
  int w,h, w2,h2;
  int dx,dy;

  w = abs(ax + ay);
  h = abs(bx + by);

  // unit major direction
  dax = sgn(ax);
  day = sgn(ay);

  // unit orthogonal direction
  dbx = sgn(bx);
  dby = sgn(by);

  dx = dax+dbx;
  dy = day+dby;

  if (h == 1) {
    if (dax==0) { return cur_idx + (dy*(y_dst-y)); }
    return cur_idx + (dx*(x_dst-x));
  }

  if (w == 1) {
    if (dbx==0) {  return cur_idx + (dy*(y_dst-y)); }
    return cur_idx + (dx*(x_dst-x));
  }

  ax2 = ax/2;
  ay2 = ay/2;
  bx2 = bx/2;
  by2 = by/2;

  w2 = abs(ax2 + ay2);
  h2 = abs(bx2 + by2);

  if ((2*w) > (3*h)) {
    if ((w2 % 2) && (w > 2)) {
      // prefer even steps
      ax2 += dax;
      ay2 += day;
    }

    if (inbounds2( x_dst, y_dst, x,y, ax2,ay2, bx,by )) {
      return gilbertxy2d_r(cur_idx, x_dst, y_dst, x, y, ax2, ay2, bx, by);
    }
    cur_idx += abs((ax2 + ay2)*(bx + by));

    return gilbertxy2d_r(cur_idx, x_dst, y_dst, x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by);
  }

  if ((h2 % 2) && (h > 2)) {
    // prefer even steps
    bx2 += dbx;
    by2 += dby;
  }

  // standard case: one step up, one long horizontal, one step down
  if (inbounds2( x_dst,y_dst, x,y, bx2,by2, ax2,ay2 )) {
    return gilbertxy2d_r(cur_idx, x_dst,y_dst, x,y, bx2,by2, ax2,ay2);
  }
  cur_idx += abs((bx2 + by2)*(ax2 + ay2));

  if (inbounds2( x_dst,y_dst, x+bx2, y+by2, ax, ay, bx-bx2, by-by2)) {
    return gilbertxy2d_r(cur_idx, x_dst,y_dst, x+bx2,y+by2, ax,ay, bx-bx2,by-by2);
  }
  cur_idx += abs((ax+ay)*((bx-bx2) + (by-by2)));

  return gilbertxy2d_r(cur_idx, x_dst,y_dst,
                 x+(ax-dax)+(bx2-dbx),
                 y+(ay-day)+(by2-dby),
                 -bx2, -by2,
                 -(ax-ax2), -(ay-ay2));
}



int gilbertd2xyz_r(int dst_idx, int cur_idx,
                   int *xres, int *yres, int *zres,
                   int ax, int ay, int az,
                   int bx, int by, int bz,
                   int cx, int cy, int cz) {
  int x,y,z;
  int _dx,_dy,_dz,_di;
  int nxt_idx;

  int w,h,d;
  int w2,h2,d2;

  int dax,day,daz,
      dbx,dby,dbz,
      dcx,dcy,dcz;
  int ax2,ay2,az2,
      bx2,by2,bz2,
      cx2,cy2,cz2;

  x = *xres;
  y = *yres;
  z = *zres;

  w = abs(ax + ay + az);
  h = abs(bx + by + bz);
  d = abs(cx + cy + cz);

  dax=sgn(ax); day=sgn(ay); daz=sgn(az);  // unit major direction "right"
  dbx=sgn(bx); dby=sgn(by); dbz=sgn(bz);  // unit ortho direction "forward"
  dcx=sgn(cx); dcy=sgn(cy); dcz=sgn(cz);  // unit ortho direction "up"

  _dx = dax+dbx+dcx;
  _dy = day+dby+dcy;
  _dz = daz+dbz+dcz;
  _di = dst_idx - cur_idx;

  // trivial row/column fills
  if ((h == 1) && (d == 1)) {
    *xres = x + dax*_di;
    *yres = y + day*_di;
    *zres = z + daz*_di;
    return 0;
  }

  if ((w == 1) && (d == 1)) {
    *xres = x + dbx*_di;
    *yres = y + dby*_di;
    *zres = z + dbz*_di;
    return 0;
  }

  if ((w == 1) && (h == 1)) {
    *xres = x + dcx*_di;
    *yres = y + dcy*_di;
    *zres = z + dcz*_di;
    return 0;
  }

  ax2=ax/2; ay2=ay/2; az2=az/2;
  bx2=bx/2; by2=by/2; bz2=bz/2;
  cx2=cx/2; cy2=cy/2; cz2=cz/2;

  w2 = abs(ax2 + ay2 + az2);
  h2 = abs(bx2 + by2 + bz2);
  d2 = abs(cx2 + cy2 + cz2);

  // prefer even steps
  if ((w2 % 2) && (w > 2)) { ax2+=dax; ay2+=day; az2+=daz; }
  if ((h2 % 2) && (h > 2)) { bx2+=dbx; by2+=dby; bz2+=dbz; }
  if ((d2 % 2) && (d > 2)) { cx2+=dcx; cy2+=dcy; cz2+=dcz; }

  // wide case, split in w only
  if (((2*w) > (3*h)) && ((2*w) > (3*d))) {
    nxt_idx = cur_idx + abs( (ax2+ay2+az2)*(bx+by+bz)*(cx+cy+cz) );
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      *xres = x;
      *yres = y;
      *zres = z;
      return gilbertd2xyz_r(dst_idx,cur_idx,
                            xres,yres,zres,
                            ax2, ay2, az2,
                            bx, by, bz,
                            cx, cy, cz);
    }
    cur_idx = nxt_idx;


    *xres = x+ax2;
    *yres = y+ay2;
    *zres = z+az2;
    return gilbertd2xyz_r(dst_idx,cur_idx,
                          xres,yres,zres,
                          ax-ax2, ay-ay2, az-az2,
                          bx, by, bz,
                          cx, cy, cz);
  }

  // do not split in d
  else if ((3*h) > (4*d)) {
    nxt_idx = cur_idx + abs( (bx2+by2+bz2)*(cx+cy+cz)*(ax2+ay2+az2) );
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      *xres = x;
      *yres = y;
      *zres = z;
      return gilbertd2xyz_r(dst_idx,cur_idx,
                            xres,yres,zres,
                            bx2, by2, bz2,
                            cx, cy, cz,
                            ax2, ay2, az2);
    }
    cur_idx = nxt_idx;

    nxt_idx = cur_idx + abs( (ax+ay+az)*((bx-bx2)+(by-by2)+(bz-bz2))*(cx+cy+cz) );
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      *xres = x+bx2;
      *yres = y+by2;
      *zres = z+bz2;
      return gilbertd2xyz_r(dst_idx,cur_idx,
                            xres,yres,zres,
                            ax, ay, az,
                            bx-bx2, by-by2, bz-bz2,
                            cx, cy, cz);
    }
    cur_idx = nxt_idx;

    *xres = x+(ax-dax)+(bx2-dbx);
    *yres = y+(ay-day)+(by2-dby);
    *zres = z+(az-daz)+(bz2-dbz);

    return gilbertd2xyz_r(dst_idx,cur_idx,
                          xres,yres,zres,
                          -bx2, -by2, -bz2,
                          cx, cy, cz,
                          -(ax-ax2), -(ay-ay2), -(az-az2));
  }

  // do not split in h
  else if ((3*d) > (4*h)) {
    nxt_idx = cur_idx + abs( (cx2+cy2+cz2)*(ax2+ay2+az2)*(bx+by+bz) );
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      *xres = x;
      *yres = y;
      *zres = z;
      return gilbertd2xyz_r(dst_idx,cur_idx,
                            xres,yres,zres,
                            cx2, cy2, cz2,
                            ax2, ay2, az2,
                            bx, by, bz);
    }
    cur_idx = nxt_idx;

    nxt_idx = cur_idx + abs( (ax+ay+az)*(bx+by+bz)*((cx-cx2)+(cy-cy2)+(cz-cz2)) );
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      *xres = x+cx2;
      *yres = y+cy2;
      *zres = z+cz2;
      return gilbertd2xyz_r(dst_idx,cur_idx,
                            xres,yres,zres,
                            ax, ay, az,
                            bx, by, bz,
                            cx-cx2, cy-cy2, cz-cz2);
    }
    cur_idx = nxt_idx;

    *xres = x+(ax-dax)+(cx2-dcx);
    *yres = y+(ay-day)+(cy2-dcy);
    *zres = z+(az-daz)+(cz2-dcz);

    return gilbertd2xyz_r(dst_idx,cur_idx,
                          xres,yres,zres,
                          -cx2, -cy2, -cz2,
                          -(ax-ax2), -(ay-ay2), -(az-az2),
                          bx, by, bz);

  }

  // regular case, split in all w/h/d
  nxt_idx = cur_idx + abs( (bx2+by2+bz2)*(cx2+cy2+cz2)*(ax2+ay2+az2) );
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    *xres = x;
    *yres = y;
    *zres = z;
    return gilbertd2xyz_r(dst_idx,cur_idx,
                          xres,yres,zres,
                          bx2, by2, bz2,
                          cx2, cy2, cz2,
                          ax2, ay2, az2);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + abs( (cx+cy+cz)*(ax2+ay2+az2)*((bx-bx2)+(by-by2)+(bz-bz2)) );
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    *xres = x+bx2;
    *yres = y+by2;
    *zres = z+bz2;
    return gilbertd2xyz_r(dst_idx,cur_idx,
                          xres,yres,zres,
                          cx, cy, cz,
                          ax2, ay2, az2,
                          bx-bx2, by-by2, bz-bz2);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + abs( (ax+ay+az)*(-bx2-by2-bz2)*(-(cx-cx2)-(cy-cy2)-(cz-cz2)) );
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    *xres = x+(bx2-dbx)+(cx-dcx);
    *yres = y+(by2-dby)+(cy-dcy);
    *zres = z+(bz2-dbz)+(cz-dcz);
    return gilbertd2xyz_r(dst_idx, cur_idx,
                          xres,yres,zres,
                          ax, ay, az,
                          -bx2, -by2, -bz2,
                          -(cx-cx2), -(cy-cy2), -(cz-cz2));
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + abs( (-cx-cy-cz)*(-(ax-ax2)-(ay-ay2)-(az-az2))*((bx-bx2)+(by-by2)+(bz-bz2)) );
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    *xres = x+(ax-dax)+bx2+(cx-dcx);
    *yres = y+(ay-day)+by2+(cy-dcy);
    *zres = z+(az-daz)+bz2+(cz-dcz);
    return gilbertd2xyz_r(dst_idx,cur_idx,
                          xres,yres,zres,
                          -cx, -cy, -cz,
                          -(ax-ax2), -(ay-ay2), -(az-az2),
                          bx-bx2, by-by2, bz-bz2);
  }
  cur_idx = nxt_idx;

  *xres = x+(ax-dax)+(bx2-dbx);
  *yres = y+(ay-day)+(by2-dby);
  *zres = z+(az-daz)+(bz2-dbz);
  return gilbertd2xyz_r(dst_idx,cur_idx,
                        xres,yres,zres,
                        -bx2, -by2, -bz2,
                        cx2, cy2, cz2,
                        -(ax-ax2), -(ay-ay2), -(az-az2));

}



int gilbertxyz2d_r(int cur_idx,
                   int x_dst, int y_dst, int z_dst,
                   int x, int y, int z,
                   int ax, int ay, int az,
                   int bx, int by, int bz,
                   int cx, int cy, int cz) {
  int w,h,d;
  int w2,h2,d2;

  int dax,day,daz,
      dbx,dby,dbz,
      dcx,dcy,dcz;

  int ax2,ay2,az2,
      bx2,by2,bz2,
      cx2,cy2,cz2;

  w = abs(ax + ay + az);
  h = abs(bx + by + bz);
  d = abs(cx + cy + cz);

  dax=sgn(ax); day=sgn(ay); daz=sgn(az); // unit major direction ("right")
  dbx=sgn(bx); dby=sgn(by); dbz=sgn(bz); // unit ortho direction ("forward")
  dcx=sgn(cx); dcy=sgn(cy); dcz=sgn(cz); // unit ortho direction ("up")

  // trivial row/column fills
  if ((h == 1) && (d == 1)) {
    return cur_idx + (dax*(x_dst-x)) + (day*(y_dst-y)) + (daz*(z_dst-z)); 
  }

  if ((w == 1) && (d == 1)) {
    return cur_idx + (dbx*(x_dst-x)) + (dby*(y_dst-y)) + (dbz*(z_dst-z)); 
  }

  if ((w == 1) && (h == 1)) {
    return cur_idx + (dcx*(x_dst-x)) + (dcy*(y_dst-y)) + (dcz*(z_dst-z)); 
  }

  ax2=ax/2; ay2=ay/2; az2=az/2;
  bx2=bx/2; by2=by/2; bz2=bz/2;
  cx2=cx/2; cy2=cy/2; cz2=cz/2;

  w2 = abs(ax2 + ay2 + az2);
  h2 = abs(bx2 + by2 + bz2);
  d2 = abs(cx2 + cy2 + cz2);

  // prefer even steps
  if ((w2 % 2) && (w > 2)) {
    ax2 += dax;
    ay2 += day;
    az2 += daz;
  }

  if ((h2 % 2) && (h > 2)) {
    bx2 += dbx;
    by2 += dby;
    bz2 += dbz;
  }

  if ((d2 % 2) && (d > 2)) {
    cx2 += dcx;
    cy2 += dcy;
    cz2 += dcz;
  }

  // wide case, split in w only
  if ((2*w > 3*h) && (2*w > 3*d)) {
    if (inbounds3(x_dst,y_dst,z_dst,
                  x,y,z,
                  ax2,ay2,az2,
                  bx,by,bz,
                  cx,cy,cz)) {
      return gilbertxyz2d_r(cur_idx,
                            x_dst,y_dst,z_dst,
                            x, y, z,
                            ax2, ay2, az2,
                            bx, by, bz,
                            cx, cy, cz);
    }
    cur_idx += abs( (ax2+ay2+az2)*(bx+by+bz)*(cx+cy+cz) );

    return gilbertxyz2d_r(cur_idx,
                          x_dst,y_dst,z_dst,
                          x+ax2, y+ay2, z+az2,
                          ax-ax2, ay-ay2, az-az2,
                          bx, by, bz,
                          cx, cy, cz);
  }

  // do not split in d
  else if ((3*h) > (4*d)) {
    if (inbounds3(x_dst,y_dst,z_dst,
                  x,y,z,
                  bx2,by2,bz2,
                  cx,cy,cz,
                  ax2,ay2,az2)) {
      return gilbertxyz2d_r(cur_idx,
                            x_dst,y_dst,z_dst,
                            x, y, z,
                            bx2, by2, bz2,
                            cx, cy, cz,
                            ax2, ay2, az2);
    }
    cur_idx += abs( (bx2+by2+bz2)*(cx+cy+cz)*(ax2+ay2+az2) );

    if (inbounds3(x_dst,y_dst,z_dst,
                  x+bx2,y+by2,z+bz2,
                  ax,ay,az,
                  bx-bx2,by-by2,bz-bz2,
                  cx,cy,cz)) {
      return gilbertxyz2d_r(cur_idx,
                            x_dst,y_dst,z_dst,
                            x+bx2, y+by2, z+bz2,
                            ax, ay, az,
                            bx-bx2, by-by2, bz-bz2,
                            cx, cy, cz);
    }
    cur_idx += abs( (ax+ay+az)*((bx-bx2)+(by-by2)+(bz-bz2))*(cx+cy+cz) );

    return gilbertxyz2d_r(cur_idx,
                          x_dst,y_dst,z_dst,
                          x+(ax-dax)+(bx2-dbx),
                          y+(ay-day)+(by2-dby),
                          z+(az-daz)+(bz2-dbz),
                          -bx2, -by2, -bz2,
                          cx, cy, cz,
                          -(ax-ax2), -(ay-ay2), -(az-az2));
  }

  // do not split in h
  else if ((3*d) > (4*h)) {
    if (inbounds3(x_dst,y_dst,z_dst,
                  x,y,z,
                  cx2,cy2,cz2,
                  ax2,ay2,az2, bx,by,bz)) {
      return gilbertxyz2d_r(cur_idx,
                            x_dst,y_dst,z_dst,
                            x, y, z,
                            cx2, cy2, cz2,
                            ax2, ay2, az2,
                            bx, by, bz);
    }
    cur_idx += abs( (cx2+cy2+cz2)*(ax2+ay2+az2)*(bx+by+bz) );

    if (inbounds3(x_dst,y_dst,z_dst,
                  x+cx2,y+cy2,z+cz2,
                  ax,ay,az, bx,by,bz,
                  cx-cx2,cy-cy2,cz-cz2)) {
      return gilbertxyz2d_r(cur_idx,
                            x_dst,y_dst,z_dst,
                            x+cx2, y+cy2, z+cz2,
                            ax, ay, az,
                            bx, by, bz,
                            cx-cx2, cy-cy2, cz-cz2);
    }
    cur_idx += abs( (ax+ay+az)*(bx+by+bz)*((cx-cx2)+(cy-cy2)+(cz-cz2)) );

    return gilbertxyz2d_r(cur_idx,
                          x_dst,y_dst,z_dst,
                          x+(ax-dax)+(cx2-dcx),
                          y+(ay-day)+(cy2-dcy),
                          z+(az-daz)+(cz2-dcz),
                          -cx2, -cy2, -cz2,
                          -(ax-ax2), -(ay-ay2), -(az-az2),
                          bx, by, bz);

  }

  // regular case, split in all w/h/d
  if (inbounds3(x_dst,y_dst,z_dst,
                x,y,z,
                bx2,by2,bz2,
                cx2,cy2,cz2,
                ax2,ay2,az2)) {
    return gilbertxyz2d_r(cur_idx,x_dst,y_dst,z_dst,
                          x, y, z,
                          bx2, by2, bz2,
                          cx2, cy2, cz2,
                          ax2, ay2, az2);
  }
  cur_idx += abs( (bx2+by2+bz2)*(cx2+cy2+cz2)*(ax2+ay2+az2) );

  if (inbounds3(x_dst,y_dst,z_dst,
                x+bx2, y+by2, z+bz2,
                cx, cy, cz,
                ax2, ay2, az2,
                bx-bx2, by-by2, bz-bz2)) {
    return gilbertxyz2d_r(cur_idx,
                          x_dst,y_dst,z_dst,
                          x+bx2, y+by2, z+bz2,
                          cx, cy, cz,
                          ax2, ay2, az2,
                          bx-bx2, by-by2, bz-bz2);
  }
  cur_idx += abs( (cx+cy+cz)*(ax2+ay2+az2)*((bx-bx2)+(by-by2)+(bz-bz2)) );

  if (inbounds3(x_dst,y_dst,z_dst,
                x+(bx2-dbx)+(cx-dcx),
                y+(by2-dby)+(cy-dcy),
                z+(bz2-dbz)+(cz-dcz),
                ax, ay, az,
                -bx2, -by2, -bz2,
                -(cx-cx2), -(cy-cy2), -(cz-cz2))) {
    return gilbertxyz2d_r(cur_idx,
                          x_dst,y_dst,z_dst,
                          x+(bx2-dbx)+(cx-dcx),
                          y+(by2-dby)+(cy-dcy),
                          z+(bz2-dbz)+(cz-dcz),
                          ax, ay, az,
                          -bx2, -by2, -bz2,
                          -(cx-cx2), -(cy-cy2), -(cz-cz2));
  }
  cur_idx += abs( (ax+ay+az)*(-bx2-by2-bz2)*(-(cx-cx2)-(cy-cy2)-(cz-cz2)) );

  if (inbounds3(x_dst,y_dst,z_dst,
                x+(ax-dax)+bx2+(cx-dcx),
                y+(ay-day)+by2+(cy-dcy),
                z+(az-daz)+bz2+(cz-dcz),
                -cx, -cy, -cz,
                -(ax-ax2), -(ay-ay2), -(az-az2),
                bx-bx2, by-by2, bz-bz2)) {
    return gilbertxyz2d_r(cur_idx,
                          x_dst,y_dst,z_dst,
                          x+(ax-dax)+bx2+(cx-dcx),
                          y+(ay-day)+by2+(cy-dcy),
                          z+(az-daz)+bz2+(cz-dcz),
                          -cx, -cy, -cz,
                          -(ax-ax2), -(ay-ay2), -(az-az2),
                          bx-bx2, by-by2, bz-bz2);
  }
  cur_idx += abs( (-cx-cy-cz)*(-(ax-ax2)-(ay-ay2)-(az-az2))*((bx-bx2)+(by-by2)+(bz-bz2)) );

  return gilbertxyz2d_r(cur_idx,
                        x_dst,y_dst,z_dst,
                        x+(ax-dax)+(bx2-dbx),
                        y+(ay-day)+(by2-dby),
                        z+(az-daz)+(bz2-dbz),
                        -bx2, -by2, -bz2,
                        cx2, cy2, cz2,
                        -(ax-ax2), -(ay-ay2), -(az-az2));
}



#define GILBERT_MAIN
#ifdef GILBERT_MAIN

#include <string.h>

int main(int argc, char **argv) {
  int w, h, d;
  int x, y, z;
  int idx;

  char buf[1024];

  w = 1;
  h = 1;
  d = 1;

  if (argc < 4) {
    printf("provide args\n");
    printf("\n");
    printf("usage:\n");
    printf("\n");
    printf("  gilbert <op> <width> <height> [depth]\n");
    printf("\n");
    printf("    op      - one of \"xy2d\",\"2dxy\",\"xyz2d\",\"d2xyz\"\n");
    printf("    depth   - default to 1 for 3D Gilbert with no depth specified\n");
    printf("\n");
    exit(-1);
  }

  strncpy(buf, argv[1], 1023);
  buf[1024]='\0';

  w = atoi(argv[2]);
  h = atoi(argv[3]);
  if (argc > 4) {
    d = atoi(argv[4]);
  }

  if ((w<=0) || (h<=0) || (d<=0)) {
    exit(-1);
  }

  if (strncmp("xy2d", buf, 1023)==0) {

    for (x=0; x<w; x++) {
      for (y=0; y<h; y++) {
        idx = gilbertxy2d( x, y, w, h );
        printf("%i %i %i\n", idx, x, y);
      }
    }

  }
  else if (strncmp("d2xy", buf, 1023)==0) {

    for (idx=0; idx<(w*h); idx++) {
      gilbertd2xy( &x, &y, idx, w, h );
      printf("%i %i\n", x, y);
    }

  }
  else if (strncmp("xyz2d", buf, 1023)==0) {

    for (x=0; x<w; x++) {
      for (y=0; y<h; y++) {
        for (z=0; z<d; z++) {
          idx = gilbertxyz2d( x,y,z, w,h,d );
          printf("%i %i %i %i\n", idx, x, y, z);
        }
      }
    }

  }

  else if (strncmp("d2xyz", buf, 1023)==0) {

    for (idx=0; idx<(w*h*d); idx++) {
      gilbertd2xyz( &x,&y,&z, idx, w,h,d );
      printf("%i %i %i\n", x, y, z);
    }

  }

  exit(0);
}

#endif
