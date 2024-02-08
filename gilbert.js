// SPDX-License-Identifier: BSD-2-Clause
// Copyright (c) 2018 Jakub Červený

"use strict";


var gilbert = {
  "xy2d": gilbertxy2d,
  "d2xy": gilbertd2xy,

  "xyz2d": gilbertxyz2d,
  "d2xyz": gilbertd2xyz,
};

function sgn(x) {
  if (x<0) { return -1; }
  if (x>0) { return  1; }
  return 0;
}

function inbounds2(p, s, a, b) {
  let d = { "x":a.x+b.x, "y":a.y+b.y };

  if (d.x < 0) {
    if ((p.x > s.x) || (p.x <= (s.x+d.x))) { return false; }
  }
  else if ((p.x < s.x) || (p.x >= (s.x+d.x))) { return false; }

  if (d.y<0) {
    if ((p.y > s.y) || (p.y <= (s.y+d.y))) { return false; }
  }
  else if ((p.y < s.y) || (p.y >= (s.y+d.y))) { return false; }

  return true;
}

function inbounds3(p, s, a, b, c) {
  let d = { "x":a.x+b.x+c.x, "y":a.y+b.y+c.y, "z":a.z+b.z+c.z };

  if (d.x < 0) {
    if ((p.x > s.x) || (p.x <= (s.x+d.x))) { return false; }
  }
  else if ((p.x < s.x) || (p.x >= (s.x+d.x))) { return false; }

  if (d.y<0) {
    if ((p.y > s.y) || (p.y <= (s.y+d.y))) { return false; }
  }
  else if ((p.y < s.y) || (p.y >= (s.y+d.y))) { return false; }

  if (d.z<0) {
    if ((p.z > s.z) || (p.z <= (s.z+d.z))) { return false; }
  }
  else if ((p.z < s.z) || (p.z >= (s.z+d.z))) { return false; }

  return true;
}

function gilbertxy2d(x,y,w,h) {
  let _q = {"x":x, "y":y};
  let _p = {"x":0, "y":0};
  let _a = {"x":0, "y":h};
  let _b = {"x":w, "y":0};

  if (w >= h) {
    _a.x = w; _a.y = 0;
    _b.x = 0; _b.y = h;
  }
  return gilbertxy2d_r(0, _q, _p, _a, _b);
}

function gilbertd2xy(idx,w,h) {
  let _p = {"x":0, "y":0};
  let _a = {"x":0, "y":h};
  let _b = {"x":w, "y":0};

  if (w >= h) {
    _a.x = w; _a.y = 0;
    _b.x = 0; _b.y = h;
  }
  return gilbertd2xy_r(idx,0,_p,_a,_b);
}

function gilbertxyz2d(x,y,z,w,h,d) {
  let _q = {"x":x, "y":y, "z":z};
  let _p = {"x":0, "y":0, "z":0};
  let _a = {"x":w, "y":0, "z": 0};
  let _b = {"x":0, "y":h, "z": 0};
  let _c = {"x":0, "y":0, "z": d};

  if ((w >= h) && (w >= d)) {
    return gilbertxyz2d_r(0, _q, _p, _a, _b, _c);
  }
  else if ((h >= w) && (h >= d)) {
    return gilbertxyz2d_r(0, _q, _p, _b, _a, _c);
  }
  return gilbertxyz2d_r(0, _q, _p, _c, _a, _b);
}

function gilbertd2xyz(idx,w,h,d) {
  let _p = {"x":0, "y":0, "z":0};
  let _a = {"x":w, "y":0, "z": 0};
  let _b = {"x":0, "y":h, "z": 0};
  let _c = {"x":0, "y":0, "z": d};
  if ((w >= h) && (w >= d)) {
    return gilbertd2xyz_r(idx, 0, _p, _a, _b, _c);
  }
  else if ((h >= w) && (h >= d)) {
    return gilbertd2xyz_r(idx, 0, _p, _b, _a, _c);
  }
  return gilbertd2xyz_r(idx, 0, _p, _c, _a, _b);
}

function gilbertd2xy_r( dst_idx,cur_idx, p, a, b) {
  let _p = {}, _a = {}, _b = {};
  let nxt_idx = -1;

  let w = Math.abs( a.x + a.y );
  let h = Math.abs( b.x + b.y );

  let da = { "x": sgn(a.x), "y": sgn(a.y) };
  let db = { "x": sgn(b.x), "y": sgn(b.y) };
  let d = { "x": da.x+db.x, "y": da.y+db.y, "i": dst_idx - cur_idx };

  if (h==1) {
    return { "x": p.x + da.x*d.i, "y": p.y + da.y*d.i }; }
  if (w==1) {
    return {"x": p.x + db.x*d.i, "y": p.y + db.y*d.i }; }

  let a2 = { "x": Math.floor(a.x/2), "y": Math.floor(a.y/2) };
  let b2 = { "x": Math.floor(b.x/2), "y": Math.floor(b.y/2) };

  let w2 = Math.abs(a2.x + a2.y);
  let h2 = Math.abs(b2.x + b2.y);


  if ((2*w) > (3*h)) {

    // prefer even steps
    if ((w2%2) && (w>2)) {
      a2.x += da.x;
      a2.y += da.y;
    }

    nxt_idx = cur_idx + Math.abs((a2.x+a2.y)*(b.x+b.y));
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return gilbertd2xy_r(dst_idx,cur_idx, p, a2, b);
    }
    cur_idx = nxt_idx;

    _p = { "x":p.x+a2.x, "y":p.y+a2.y };
    _a = { "x":a.x-a2.x, "y":a.y-a2.y };

    return gilbertd2xy_r(dst_idx,cur_idx, _p, _a, b);
  }

  // prefer event steps
  if ((h2%2) && (h>2)) {
    b2.x += db.x;
    b2.y += db.y;
  }

  // standard case: one step up, on long horizontal, one step down
  nxt_idx = cur_idx + Math.abs((b2.x+b2.y)*(a2.x+a2.y));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return gilbertd2xy_r(dst_idx, cur_idx, p, b2, a2);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + Math.abs((a.x+a.y)*((b.x-b2.x) + (b.y-b2.y)));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    _p = { "x":p.x+b2.x, "y":p.y+b2.y };
    _b = { "x":b.x-b2.x, "y":b.y-b2.y };
    return gilbertd2xy_r(dst_idx, cur_idx, _p, a, _b);
  }
  cur_idx = nxt_idx;

  _p = {
    "x":p.x+(a.x-da.x)+(b2.x-db.x),
    "y":p.y+(a.y-da.y)+(b2.y-db.y)
  };
  _a = { "x": -b2.x, "y": -b2.y };
  _b = { "x": -(a.x-a2.x), "y": -(a.y-a2.y) };

  return gilbertd2xy_r(dst_idx, cur_idx, _p, _a, _b);
}

function gilbertxy2d_r(idx, q, p, a, b) {
  let _p = {}, _a = {}, _b = {};

  let w = Math.abs(a.x+a.y);
  let h = Math.abs(b.x+b.y);

  let da = { "x": sgn(a.x), "y": sgn(a.y) };
  let db = { "x": sgn(b.x), "y": sgn(b.y) };
  let d = {"x": da.x+db.x, "y": da.y+db.y };

  if (h==1) {
    return idx + (da.x*(q.x-p.x)) + (da.y*(q.y-p.y));
  }
  if (w==1) {
    return idx + (db.x*(q.x-p.x)) + (db.y*(q.y-p.y));
  }

  let a2 = { "x":Math.floor(a.x/2), "y":Math.floor(a.y/2) };
  let b2 = { "x":Math.floor(b.x/2), "y":Math.floor(b.y/2) };

  let w2 = Math.abs(a2.x+a2.y);
  let h2 = Math.abs(b2.x+b2.y);

  if ((2*w) > (3*h)) {
    if ((w2%2) && (w>2)) {
      a2.x += da.x;
      a2.y += da.y;
    }

    if (inbounds2( q, p, a2, b )) {
      return gilbertxy2d_r(idx, q, p, a2, b);
    }
    idx += Math.abs((a2.x+a2.y)*(b.x+b.y));

    _p = { "x": p.x+a2.x, "y": p.y+a2.y };
    _a = { "x": a.x-a2.x, "y": a.y-a2.y };
    return gilbertxy2d_r(idx, q, _p, _a, b);
  }

  if ((h2%2) && (h>2)) {
    b2.x += db.x;
    b2.y += db.y;
  }

  if (inbounds2( q, p, b2, a2 )) {
    return gilbertxy2d_r(idx, q, p, b2, a2);
  }
  idx += Math.abs((b2.x+b2.y)*(a2.x+a2.y));

  _p = { "x": p.x+b2.x, "y": p.y+b2.y };
  _b = { "x": b.x-b2.x, "y": b.y-b2.y };
  if (inbounds2( q, _p, a, _b )) {
    return gilbertxy2d_r(idx, q, _p, a, _b);
  }
  idx += Math.abs((a.x+a.y)*((b.x-b2.x) + (b.y-b2.y)));

  _p = {
    "x" : p.x+(a.x-da.x)+(b2.x-db.x),
    "y" : p.y+(a.y-da.y)+(b2.y-db.y)
  };
  _a = { "x": -b2.x, "y": -b2.y };
  _b = { "x": -(a.x-a2.x), "y": -(a.y-a2.y) };
  return gilbertxy2d_r(idx, q, _p, _a, _b);

}


function gilbertxyz2d_r(cur_idx, q, p, a, b, c) {
  let _p = {}, _a = {}, _b = {}, _c = {};

  let w = Math.abs(a.x+a.y+a.z);
  let h = Math.abs(b.x+b.y+b.z);
  let d = Math.abs(c.x+c.y+c.z);

  let da = { "x": sgn(a.x), "y": sgn(a.y), "z": sgn(a.z) };
  let db = { "x": sgn(b.x), "y": sgn(b.y), "z": sgn(b.z) };
  let dc = { "x": sgn(c.x), "y": sgn(c.y), "z": sgn(c.z) };

  // trivial row/column fills
  if ((h==1) && (d==1)) {
    return cur_idx + (da.x*(q.x-p.x)) + (da.y*(q.y-p.y)) + (da.z*(q.z-p.z));
  }
  else if ((w==1) && (d==1)) {
    return cur_idx + (db.x*(q.x-p.x)) + (db.y*(q.y-p.y)) + (db.z*(q.z-p.z));
  }
  else if ((w==1) && (h==1)) {
    return cur_idx + (dc.x*(q.x-p.x)) + (dc.y*(q.y-p.y)) + (dc.z*(q.z-p.z));
  }

  let a2 = { "x":Math.floor(a.x/2), "y":Math.floor(a.y/2), "z":Math.floor(a.z/2) };
  let b2 = { "x":Math.floor(b.x/2), "y":Math.floor(b.y/2), "z":Math.floor(b.z/2) };
  let c2 = { "x":Math.floor(c.x/2), "y":Math.floor(c.y/2), "z":Math.floor(c.z/2) };

  let w2 = Math.abs(a2.x+a2.y+a2.z);
  let h2 = Math.abs(b2.x+b2.y+b2.z);
  let d2 = Math.abs(c2.x+c2.y+c2.z);

  // prefer even steps
  if ((w2%2) && (w>2)) {
    a2.x += da.x;
    a2.y += da.y;
    a2.z += da.z;
  }

  if ((h2%2) && (h>2)) {
    b2.x += db.x;
    b2.y += db.y;
    b2.z += db.z;
  }

  if ((d2%2) && (d>2)) {
    c2.x += dc.x;
    c2.y += dc.y;
    c2.z += dc.z;
  }

  // wide case, split in w only
  if ( ((2*w) > (3*h)) && ((2*w) > (3*d)) ) {
    if (inbounds3(q,p,a2,b,c)) {
      return gilbertxyz2d_r(cur_idx, q, p, a2, b, c);
    }
    cur_idx += Math.abs( (a2.x+a2.y+a2.z)*(b.x+b.y+b.z)*(c.x+c.y+c.z) );

    _p = { "x":p.x+a2.x, "y":p.y+a2.y, "z":p.z+a2.z };
    _a = { "x":a.x-a2.x, "y":a.y-a2.y, "z":a.z-a2.z };
    return gilbertxyz2d_r(cur_idx, q, _p, _a, b, c);
  }

  else if ((3*h) > (4*d)) {

    if (inbounds3(q,p,b2,c,a2)) {
      return gilbertxyz2d_r(cur_idx,q,p,b2,c,a2);
    }
    cur_idx += Math.abs( (b2.x+b2.y+b2.z)*(c.x+c.y+c.z)*(a2.x+a2.y+a2.z) );

    _p = { "x":p.x+b2.x, "y":p.y+b2.y, "z":p.z+b2.z };
    _b = { "x":b.x-b2.x, "y":b.y-b2.y, "z":b.z-b2.z };
    if (inbounds3(q,_p,a,_b,c)) {
      return gilbertxyz2d_r(cur_idx,q, _p, a, _b, c);
    }
    cur_idx += Math.abs( (a.x+a.y+a.z)*((b.x-b2.x)+(b.y-b2.y)+(b.z-b2.z))*(c.x+c.y+c.z) );

    _p = {
      "x": p.x+(a.x-da.x)+(b2.x-db.x),
      "y": p.y+(a.y-da.y)+(b2.y-db.y),
      "z": p.z+(a.z-da.z)+(b2.z-db.z)
    };
    _a = { "x":-b2.x, "y":-b2.y, "z":-b2.z };
    _c = { "x":-(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
    return gilbertxyz2d_r(cur_idx, q, _p, _a, c, _c);

  }

  else if ((3*d) > (4*h)) {

    if (inbounds3(q,p,c2,a2,b)) {
      return gilbertxyz2d_r(cur_idx,q,p,c2,a2,b);
    }
    cur_idx += Math.abs( (c2.x+c2.y+c2.z)*(a2.x+a2.y+a2.z)*(b.x+b.y+b.z) );

    _p = { "x": p.x+c2.x, "y": p.y+c2.y, "z": p.z+c2.z };
    _c = { "x": c.x-c2.x, "y": c.y-c2.y, "z": c.z-c2.z };
    if (inbounds3(q, _p, a, b, _c)) {
      return gilbertxyz2d_r(cur_idx, q, _p, a, b, _c);
    }
    cur_idx += Math.abs( (a.x+a.y+a.z)*(b.x+b.y+b.z)*((c.x-c2.x)+(c.y-c2.y)+(c.z-c2.z)) );

    _p = {
      "x": p.x+(a.x-da.x)+(c2.x-dc.x),
      "y": p.y+(a.y-da.y)+(c2.y-dc.y),
      "z": p.z+(a.z-da.z)+(c2.z-dc.z)
    }
    _a = { "x": -c2.x, "y": -c2.y, "z": -c2.z };
    _b = { "x": -(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
    return gilbertxyz2d_r(cur_idx, q, _p, _a, _b, b);

  }

  // regular case, split in all w/h/d
  if (inbounds3(q,p,b2,c2,a2)) {
    return gilbertxyz2d_r(cur_idx,q,p,b2,c2,a2);
  }
  cur_idx += Math.abs( (b2.x+b2.y+b2.z)*(c2.x+c2.y+c2.z)*(a2.x+a2.y+a2.z) );

  _p = { "x": p.x+b2.x, "y": p.y+b2.y, "z": p.z+b2.z };
  _c = { "x": b.x-b2.x, "y": b.y-b2.y, "z": b.z-b2.z };
  if (inbounds3(q, _p, c, a2, _c)) {
    return gilbertxyz2d_r(cur_idx, q, _p, c, a2, _c);
  }
  cur_idx += Math.abs( (c.x+c.y+c.z)*(a2.x+a2.y+a2.z)*((b.x-b2.x)+(b.y-b2.y)+(b.z-b2.z)) );

  _p = {
    "x" : p.x+(b2.x-db.x)+(c.x-dc.x),
    "y" : p.y+(b2.y-db.y)+(c.y-dc.y),
    "z" : p.z+(b2.z-db.z)+(c.z-dc.z)
  };
  _b = { "x": -b2.x, "y": -b2.y, "z": -b2.z };
  _c = { "x": -(c.x-c2.x), "y": -(c.y-c2.y), "z": -(c.z-c2.z) };
  if (inbounds3(q, _p, a, _b, _c)) {
    return gilbertxyz2d_r(cur_idx, q, _p, a, _b, _c);
  }
  cur_idx += Math.abs( (a.x+a.y+a.z)*(-b2.x-b2.y-b2.z)*(-(c.x-c2.x)-(c.y-c2.y)-(c.z-c2.z)) );

  _p = {
    "x": p.x+(a.x-da.x)+b2.x+(c.x-dc.x),
    "y": p.y+(a.y-da.y)+b2.y+(c.y-dc.y),
    "z": p.z+(a.z-da.z)+b2.z+(c.z-dc.z)
  };
  _a = { "x": -c.x, "y": -c.y, "z": -c.z };
  _b = { "x": -(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
  _c = { "x": b.x-b2.x, "y": b.y-b2.y, "z": b.z-b2.z };
  if (inbounds3(q, _p, _a, _b, _c)) {
    return gilbertxyz2d_r(cur_idx, q, _p, _a, _b, _c);
  }
  cur_idx += Math.abs( (-c.x-c.y-c.z)*(-(a.x-a2.x)-(a.y-a2.y)-(a.z-a2.z))*((b.x-b2.x)+(b.y-b2.y)+(b.z-b2.z)) );

  _p = {
    "x": p.x+(a.x-da.x)+(b2.x-db.x),
    "y": p.y+(a.y-da.y)+(b2.y-db.y),
    "z": p.z+(a.z-da.z)+(b2.z-db.z)
  };
  _a = { "x": -b2.x, "y": -b2.y, "z": -b2.z };
  _c = { "x": -(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
  return gilbertxyz2d_r(cur_idx, q, _p, _a, c2, _c);

}

function gilbertd2xyz_r(dst_idx, cur_idx, p, a, b, c) {
  let _p = {}, _a = {}, _b = {}, _c = {};
  let nxt_idx=-1;

  let w = Math.abs(a.x+a.y+a.z);
  let h = Math.abs(b.x+b.y+b.z);
  let d = Math.abs(c.x+c.y+c.z);

  let da = { "x": sgn(a.x), "y": sgn(a.y), "z": sgn(a.z) };
  let db = { "x": sgn(b.x), "y": sgn(b.y), "z": sgn(b.z) };
  let dc = { "x": sgn(c.x), "y": sgn(c.y), "z": sgn(c.z) };
  let di = dst_idx - cur_idx;

  // trivial row/column fills
  if ((h==1) && (d==1)) {
    return { "x": p.x + da.x*di, "y": p.y + da.y*di, "z": p.z + da.z*di };
  }
  else if ((w==1) && (d==1)) {
    return { "x": p.x + db.x*di, "y": p.y + db.y*di, "z": p.z + db.z*di };
  }
  else if ((w==1) && (h==1)) {
    return { "x": p.x + dc.x*di, "y": p.y + dc.y*di, "z": p.z + dc.z*di };
  }

  let a2 = { "x":Math.floor(a.x/2), "y":Math.floor(a.y/2), "z":Math.floor(a.z/2) };
  let b2 = { "x":Math.floor(b.x/2), "y":Math.floor(b.y/2), "z":Math.floor(b.z/2) };
  let c2 = { "x":Math.floor(c.x/2), "y":Math.floor(c.y/2), "z":Math.floor(c.z/2) };

  let w2 = Math.abs(a2.x+a2.y+a2.z);
  let h2 = Math.abs(b2.x+b2.y+b2.z);
  let d2 = Math.abs(c2.x+c2.y+c2.z);

  // prefer even steps
  if ((w2%2) && (w>2)) {
    a2.x += da.x;
    a2.y += da.y;
    a2.z += da.z;
  }

  if ((h2%2) && (h>2)) {
    b2.x += db.x;
    b2.y += db.y;
    b2.z += db.z;
  }

  if ((d2%2) && (d>2)) {
    c2.x += dc.x;
    c2.y += dc.y;
    c2.z += dc.z;
  }

  // wide case, split in w only
  if ( ((2*w) > (3*h)) && ((2*w) > (3*d)) ) {
    nxt_idx = cur_idx + Math.abs( (a2.x+a2.y+a2.z)*(b.x+b.y+b.z)*(c.x+c.y+c.z) );
    if ((cur_idx <= nxt_idx) && (dst_idx < nxt_idx)) {
      return gilbertd2xyz_r(dst_idx, cur_idx, p, a2, b, c);
    }
    cur_idx = nxt_idx;

    _p = { "x":p.x+a2.x, "y":p.y+a2.y, "z":p.z+a2.z };
    _a = { "x":a.x-a2.x, "y":a.y-a2.y, "z":a.z-a2.z };
    return gilbertd2xyz_r(dst_idx, cur_idx, _p, _a, b, c);
  }

  else if ((3*h) > (4*d)) {

    nxt_idx = cur_idx + Math.abs( (b2.x+b2.y+b2.z)*(c.x+c.y+c.z)*(a2.x+a2.y+a2.z) );
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return gilbertd2xyz_r(dst_idx,cur_idx,p,b2,c,a2);
    }
    cur_idx = nxt_idx;

    nxt_idx = cur_idx + Math.abs( (a.x+a.y+a.z)*((b.x-b2.x)+(b.y-b2.y)+(b.z-b2.z))*(c.x+c.y+c.z) );
    _p = { "x":p.x+b2.x, "y":p.y+b2.y, "z":p.z+b2.z };
    _b = { "x":b.x-b2.x, "y":b.y-b2.y, "z":b.z-b2.z };
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return gilbertd2xyz_r(dst_idx,cur_idx, _p, a, _b, c);
    }
    cur_idx = nxt_idx;

    _p = {
      "x": p.x+(a.x-da.x)+(b2.x-db.x),
      "y": p.y+(a.y-da.y)+(b2.y-db.y),
      "z": p.z+(a.z-da.z)+(b2.z-db.z)
    };
    _a = { "x":-b2.x, "y":-b2.y, "z":-b2.z };
    _c = { "x":-(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
    return gilbertd2xyz_r(dst_idx, cur_idx, _p, _a, c, _c);

  }

  else if ((3*d) > (4*h)) {

    nxt_idx = cur_idx + Math.abs( (c2.x+c2.y+c2.z)*(a2.x+a2.y+a2.z)*(b.x+b.y+b.z) );
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return gilbertd2xyz_r(dst_idx,cur_idx,p,c2,a2,b);
    }
    cur_idx = nxt_idx;

    nxt_idx = cur_idx + Math.abs( (a.x+a.y+a.z)*(b.x+b.y+b.z)*((c.x-c2.x)+(c.y-c2.y)+(c.z-c2.z)) );
    _p = { "x": p.x+c2.x, "y": p.y+c2.y, "z": p.z+c2.z };
    _c = { "x": c.x-c2.x, "y": c.y-c2.y, "z": c.z-c2.z };
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return gilbertd2xyz_r(dst_idx, cur_idx, _p, a, b, _c);
    }
    cur_idx = nxt_idx;

    _p = {
      "x": p.x+(a.x-da.x)+(c2.x-dc.x),
      "y": p.y+(a.y-da.y)+(c2.y-dc.y),
      "z": p.z+(a.z-da.z)+(c2.z-dc.z)
    }
    _a = { "x": -c2.x, "y": -c2.y, "z": -c2.z };
    _b = { "x": -(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
    return gilbertd2xyz_r(dst_idx, cur_idx, _p, _a, _b, b);

  }

  // regular case, split in all w/h/d
  nxt_idx = cur_idx + Math.abs( (b2.x+b2.y+b2.z)*(c2.x+c2.y+c2.z)*(a2.x+a2.y+a2.z) );
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return gilbertd2xyz_r(dst_idx,cur_idx,p,b2,c2,a2);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + Math.abs( (c.x+c.y+c.z)*(a2.x+a2.y+a2.z)*((b.x-b2.x)+(b.y-b2.y)+(b.z-b2.z)) );
  _p = { "x": p.x+b2.x, "y": p.y+b2.y, "z": p.z+b2.z };
  _c = { "x": b.x-b2.x, "y": b.y-b2.y, "z": b.z-b2.z };
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return gilbertd2xyz_r(dst_idx, cur_idx, _p, c, a2, _c);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + Math.abs( (a.x+a.y+a.z)*(-b2.x-b2.y-b2.z)*(-(c.x-c2.x)-(c.y-c2.y)-(c.z-c2.z)) );
  _p = {
    "x" : p.x+(b2.x-db.x)+(c.x-dc.x),
    "y" : p.y+(b2.y-db.y)+(c.y-dc.y),
    "z" : p.z+(b2.z-db.z)+(c.z-dc.z)
  };
  _b = { "x": -b2.x, "y": -b2.y, "z": -b2.z };
  _c = { "x": -(c.x-c2.x), "y": -(c.y-c2.y), "z": -(c.z-c2.z) };
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return gilbertd2xyz_r(dst_idx, cur_idx, _p, a, _b, _c);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + Math.abs( (-c.x-c.y-c.z)*(-(a.x-a2.x)-(a.y-a2.y)-(a.z-a2.z))*((b.x-b2.x)+(b.y-b2.y)+(b.z-b2.z)) );
  _p = {
    "x": p.x+(a.x-da.x)+b2.x+(c.x-dc.x),
    "y": p.y+(a.y-da.y)+b2.y+(c.y-dc.y),
    "z": p.z+(a.z-da.z)+b2.z+(c.z-dc.z)
  };
  _a = { "x": -c.x, "y": -c.y, "z": -c.z };
  _b = { "x": -(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
  _c = { "x": b.x-b2.x, "y": b.y-b2.y, "z": b.z-b2.z };
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return gilbertd2xyz_r(dst_idx, cur_idx, _p, _a, _b, _c);
  }
  cur_idx = nxt_idx;

  _p = {
    "x": p.x+(a.x-da.x)+(b2.x-db.x),
    "y": p.y+(a.y-da.y)+(b2.y-db.y),
    "z": p.z+(a.z-da.z)+(b2.z-db.z)
  };
  _a = { "x": -b2.x, "y": -b2.y, "z": -b2.z };
  _c = { "x": -(a.x-a2.x), "y": -(a.y-a2.y), "z": -(a.z-a2.z) };
  return gilbertd2xyz_r(dst_idx, cur_idx, _p, _a, c2, _c);


}

if (typeof module !== "undefined") {

  module.exports["d2xy"] = gilbertxy2d;
  module.exports["xy2d"] = gilbertd2xy;

  module.exports["d2xyz"] = gilbertxyz2d;
  module.exports["xyz2d"] = gilbertd2xyz;

  module.exports["main"] = _main;

  function _main(argv) {

    if (argv.length < 4) {
      console.log("provide args");
      process.exit(-1);
    }

    let op = argv[1];
    let w = parseInt(argv[2]);
    let h = parseInt(argv[3]);
    let d = 1;
    if (argv.length > 4) {
      d = parseInt(argv[4]);
    }

    if (op == "xy2d") {
      for (let x=0; x<w; x++) {
        for (let y=0; y<h; y++) {
          let idx = gilbertxy2d(x,y,w,h);
          console.log(idx, x, y);
        }
      }
    }

    else if (op == "d2xy") {
      let n = w*h;
      for (let idx=0; idx<n; idx++) {
        let xy = gilbertd2xy(idx,w,h);
        console.log(xy.x,xy.y);
      }
    }

    else if (op == "xyz2d") {

      for (let x=0; x<w; x++) {
        for (let y=0; y<h; y++) {
          for (let z=0; z<d; z++) {
            let idx = gilbertxyz2d(x,y,z,w,h,d);
            console.log(idx,x,y,z);
          }
        }
      }

    }

    else if (op == "d2xyz") {
      let n = w*h*d;
      for (let idx=0; idx<n; idx++) {
        let xyz = gilbertd2xyz(idx,w,h,d);
        console.log(xyz.x,xyz.y,xyz.z);
      }
    }

  }

  //_main( process.argv.slice(1) );

}


