

function sgn(x) {
  if (x<0) { return -1; }
  if (x>0) { return  1; }
  return 0;
}

function inbounds2(p, s, a, b) {
  d = { "x":a.x+b.x, "y":a.y+b.y };

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

function gilbertxy2d(x,y,w,h) {
  let _q = {"x":x, "y":y};
  let _p = {"x":0, "y":0};
  let _a = {"x":w, "y":0};
  let _b = {"x":0, "y":h};
  return gilbertxy2d_r(0, _q, _p, _a, _b);
}

function gilbertd2xy(idx,w,h) {
  let _p = {"x":0, "y":0};
  let _a = {"x":w, "y":0};
  let _b = {"x":0, "y":h};
  return gilbertd2xy_r(idx,0,_p,_a,_b);
}

function gilbertd2xy_r( dst_idx,cur_idx, p, a, b) {
  let w = Math.abs( a.x + a.y );
  let h = Math.abs( b.x + b.y );

  let da = { "x": sgn(a.x), "y": sgn(a.y) };
  let db = { "x": sgn(b.x), "y": sgn(b.y) };
  let d = { "x": da.x+db.x, "y": da.y+db.y, "i": dst_idx - cur_idx };

  if (h==1) { return { "x": p.x + da.x*d.i, "y": p.y + da.y*d.i }; }
  if (w==1) { return { "x": p.x + db.x*d.i, "y": p.y + db.y*d.i }; }

  let a2 = { "x": Math.floor(a.x/2), "y": Math.floor(a.y/2) };
  let b2 = { "x": Math.floor(b.x/2), "y": Math.floor(b.y/2) };

  let w2 = Math.abs(a2.x + a2.y);
  let h2 = Math.abs(b2.x + b2.y);

  if (2*w > 3*h) {

    // prefer even steps
    if ((w2%2) && (w>2)) {
      a2.x += da.x;
      a2.y += da.y;
    }

    let nxt_idx = cur_idx + Math.abs((a2.x+a2.y)*(b.x+b.y));
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return gilbertd2xy_r(dst_idx,cur_idx, p, a2, b);
    }
    cur_idx = nxt_idx;

    let _p = { "x":p.x+a2.x, "y":p.y+a2.y };
    let _a = { "x":a.x-a2.x, "y":a.y-a2.y };
    return gilbertd2xy_r(dst_idx,cur_idx, _p, _a, b);
  }

  // prefer event steps
  if ((h2%2) && (h>2)) {
    b2.x += db.x;
    b2.y += db.y;
  }

  // standard case: one step up, on long horizontal, one step down
  let nxt_idx = cur_idx + Math.abs((b2.x+b2.y)*(a2.x+a2.y));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return gilbertd2xy_r(dst_idx, cur_idx, p, b2, a2);
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + Math.abs((a.x+a.y)*((b.x-b2.x) + (b.y-b2.y)));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    let _p = { "x":p.x+b2.x, "y":p.y+b2.y };
    let _b = { "x":b.x-b2.x, "y":b.y-b2.y };
    return gilbertd2xy_r(dst_idx, cur_idx, _p, a, _b);
  }
  cur_idx = nxt_idx;

  let _p = {
    "x":p.x+(a.x-da.x)+(b2.x-db.x),
    "y":p.y+(a.y-da.y)+(b2.y-db.y)
  };
  let _a = { "x": -b2.x, "y": -b2.y };
  let _b = { "x": -(a.x-a2.x), "y": -(a.y-a2.y) };
  return gilbertd2xy_r(dst_idx, cur_idx, _p, _a, _b);
}

function gilbertxy2d_r(idx, q, p, a, b) {

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

    let _p = { "x": p.x+a2.x, "y": p.y+a2.y };
    let _a = { "x": a.x-a2.x, "y": a.y-a2.y };
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

  let _p = { "x": p.x+b2.x, "y": p.y+b2.y };
  let _b = { "x": b.x-b2.x, "y": b.y-b2.y };
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

if (typeof module !== "undefined") {

  module.exports["d2xy"] = gilbertxy2d;
  module.exports["xy2d"] = gilbertd2xy;

  /*
  let w= 15;
  let h = 14;

  for (let x=0; x<w; x++) {
    for (let y=0; y<h; y++) {
      let idx = gilbertxy2d(x,y,w,h);
      console.log(idx,x,y);
    }
  }

  process.exit();

  for (let idx=0; idx<w*h; idx++) {
    let xy = gilbertd2xy(idx,w,h);
    console.log(xy.x, xy.y);
  }
  */
}


