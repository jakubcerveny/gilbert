// To the extent possible under law, the person who associated CC0 with
// this project has waived all copyright and related or neighboring rights
// to this project.
// 
// You should have received a copy of the CC0 legalcode along with this
// work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//

// This is a accompanying reference implementation for the
// Gilbert2D and Gilbert3D functions referenced in the
// paper https://github.com/jakubcerveny/gilbert-paper
//
// The code here's main focus is legibility, not execution
// speed or memory usage.
// There are both asynchronous functions that closely match
// the psuedo code from the paper as well as synchronous
// functions with their corresponding d2xyz and xyz2d
// functions.
//
// Running from the command line will list options.
//

var DEBUG = 0;

function _dprint() {
  let _debug = ((typeof DEBUG === "undefined") ? false : DEBUG );
  if (_debug) {
    console.log.apply(this,Array.prototype.slice.call(arguments));
  }
}

var GILBERT_ADAPT_METHOD = {
  "HARMONY": 0,
  "HAMILTONIAN" : 1,
  "AXIS" : 2
};

function Gilbert3DAdapt_d2xyz(idx, w,h,d, adapt_method) {
  adapt_method = ((typeof adapt_method === "undefined") ? 0 : adapt_method);
  let w0 = (w%2);
  let h0 = (h%2);
  let d0 = (d%2);

  let p0 = [0,0,0];

  // prioritize harmonious split
  //
  if (adapt_method == GILBERT_ADAPT_METHOD["HARMONY"]) {
    if      ((w >= h) && (w >= d)) { return Gilbert3D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,d]); }
    else if ((h >= w) && (h >= d)) { return Gilbert3D_d2xyz(idx, 0, p0, [0,h,0], [w,0,0], [0,0,d]); }
    return Gilbert3D_d2xyz(idx, 0, p0, [0,0,d], [w,0,0], [0,h,0]);
  }

  // prioritize no notch hamiltonian path
  //
  else if (adapt_method == GILBERT_ADAPT_METHOD["HAMILTONIAN"]) {

    if (w0 == 0) { return Gilbert3D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,d]); }
    if (h0 == 0) { return Gilbert3D_d2xyz(idx, 0, p0, [0,h,0], [w,0,0], [0,0,d]); }
    if (d0 == 0) { return Gilbert3D_d2xyz(idx, 0, p0, [0,0,d], [w,0,0], [0,h,0]); }
    return Gilbert3D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,d]);

  }

  // AXIS (explicit axis order)
  //
  return Gilbert3D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,d]);

}

function Gilbert3DAdapt_xyz2d(q, w,h,d, adapt_method) {
  adapt_method = ((typeof adapt_method === "undefined") ? 0 : adapt_method);
  let w0 = (w%2);
  let h0 = (h%2);
  let d0 = (d%2);

  // prioritize harmonious split
  //
  if (adapt_method == GILBERT_ADAPT_METHOD["HARMONY"]) {
    if      ((w >= h) && (w >= d)) { return Gilbert3D_xyz2d(0, q, [0,0,0], [w,0,0], [0,h,0], [0,0,d]); }
    else if ((h >= w) && (h >= d)) { return Gilbert3D_xyz2d(0, q, [0,0,0], [0,h,0], [w,0,0], [0,0,d]); }
    return Gilbert3D_xyz2d(0, q, [0,0,0], [0,0,d], [w,0,0], [0,h,0]);
  }

  // prioritize no notch hamiltonian path
  //
  else if (adapt_method == GILBERT_ADAPT_METHOD["HAMILTONIAN"]) {

    if (w0 == 0) { return Gilbert3D_xyz2d(0, q, [0,0,0], [w,0,0], [0,h,0], [0,0,d]); }
    if (h0 == 0) { return Gilbert3D_xyz2d(0, q, [0,0,0], [0,h,0], [w,0,0], [0,0,d]); }
    if (d0 == 0) { return Gilbert3D_xyz2d(0, q, [0,0,0], [0,0,d], [w,0,0], [0,h,0]); }
    return Gilbert3D_xyz2d(0, q, [0,0,0], [w,0,0], [0,h,0], [0,0,d]);

  }

  // AXIS (explicit axis order)
  //
  return Gilbert3D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,d]);

}


function Gilbert2DAdapt_d2xy(idx, w,h, adapt_method) {
  adapt_method = ((typeof adapt_method === "undefined") ? 0 : adapt_method);
  let w0 = (w%2);
  let h0 = (h%2);

  let p0 = [0,0,0];

  // prioritize harmonious split
  //
  if (adapt_method == GILBERT_ADAPT_METHOD["HARMONY"]) {

    if (w >= h) {
      let xyz = Gilbert2D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,1]);
      return [xyz[0], xyz[1]];
    }

    let xyz = Gilbert2D_d2xyz(idx, 0, p0, [0,h,0], [w,0,0], [0,0,1]);
    return [xyz[0], xyz[1]];

  }

  // prioritize no notch hamiltonian path
  //
  else if (adapt_method == GILBERT_ADAPT_METHOD["HAMILTONIAN"]) {
    if (w0 == 0) {
      let xyz = Gilbert2D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,1]);
      return [xyz[0], xyz[1]];
    }

    if (h0 == 0) {
      let xyz = Gilbert2D_d2xyz(idx, 0, p0, [0,h,0], [w,0,0], [0,0,1]);
      return [xyz[0], xyz[1]];
    }

    let xyz = Gilbert2D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,1]);
    return [xyz[0], xyz[1]];
  }

  // AXIS (explicit axis order)
  //
  let xyz = Gilbert2D_d2xyz(idx, 0, p0, [w,0,0], [0,h,0], [0,0,1]);
  return [xyz[0], xyz[1]];

}

function Gilbert2DAdapt_xy2d(_q, w,h, adapt_method) {
  adapt_method = ((typeof adapt_method === "undefined") ? 0 : adapt_method);
  let w0 = (w%2);
  let h0 = (h%2);

  let q = [_q[0],_q[1],0];
  let p0 = [0,0,0];

  // prioritize harmonious split
  //
  if (adapt_method == GILBERT_ADAPT_METHOD["HARMONY"]) {
    if (w >= h) { return Gilbert2D_xyz2d(0, q, p0, [w,0,0], [0,h,0], [0,0,1]); }
    return Gilbert2D_xyz2d(0, q, p0, [0,h,0], [w,0,0], [0,0,1]);
  }

  // prioritize no notch hamiltonian path
  //
  else if (adapt_method == GILBERT_ADAPT_METHOD["HAMILTONIAN"]) {
    if (w0 == 0) { return Gilbert2D_xyz2d(0, q, p0, [w,0,0], [0,h,0], [0,0,1]); }
    if (h0 == 0) { return Gilbert2D_xyz2d(0, q, p0, [0,h,0], [w,0,0], [0,0,1]); }
    return Gilbert2D_xyz2d(0, q, p0, [w,0,0], [0,h,0], [0,0,1]);
  }

  // AXIS (explicit axis order)
  //
  return Gilbert2D_xyz2d(0, q, p0, [w,0,0], [0,h,0], [0,0,1]);

}

// Description:
//
// d2e    : divide by 2 but force even (by adding one if need be)
// d2u    : divide by 2 but force odd (by adding one if need be)
// dqe    : divide by q but force even (by adding one if need be)
// dqu    : divide by q but force odd (by adding one if need be)
// _divq  : integer divide vector entries by q ( sgn(v) floor(|v|/q) )
// _div2  : _divq(2)
// _sgn   : return -1,0,1 if value is <0, 0 or >0 respectively
// _neg   : negate vector (multiply each entry by -1)
// _add   : add two vectors
// _abs   : absolute value all vector entries
// _delta : replace each element in vector by _sgn(.)
// _print : print vector (debugging)
// _clone : make a clone of the vector
//
// all vector operations can be done with 2d or 3d vectors (only)
//

// floor version
//
//function _d2e(v) {
//  let v2 = Math.floor(v/2);
//  if (v==0) { return 0; }
//  if ((v2%2)==0) { return v2; }
//  return v2+1;
//}

// round to zero version
//
function d2e(_v) {
  let v = Math.abs(_v);
  let m = ((_v<0)? -1 : 1);
  let v2 = Math.floor(v/2);
  if (v==0) { return 0; }
  if ((v2%2)==0) { return m*v2; }
  return m*(v2+1);
}

// floor version
//
//function d2u_floor(v) {
//  let v2 = Math.floor(v/2);
//  if (v==0) { return 0; }
//  if ((v2%2)==1) { return v2; }
//  return v2+1;
//}

// round to zero version
//
function d2u(_v) {
  let v = Math.abs(_v);
  let m = ((_v<0)? -1 : 1);
  let v2 = Math.floor(v/2);
  if (v==0) { return 0; }
  if ((v2%2)==1) { return m*v2; }
  return m*(v2+1);
}

// floor version
//
//function dqe_floor(v,q) {
//  let vq = Math.floor(v/q);
//  if ((vq%2)==0) { return vq; }
//  return vq+1;
//}

// round to zero version
//
function dqe(_v,q) {
  let v = Math.abs(_v);
  let m = ((_v<0)? -1 : 1);
  let vq = Math.floor(v/q);
  if ((vq%2)==0) { return m*vq; }
  return m*(vq+1);
}

// floor version
//
//function dqu_floor(v,q) {
//  let vq = Math.floor(v/q);
//  if (v==0) { return 0; }
//  if ((vq%2)==1) { return vq; }
//  return vq+1;
//}


function dqu(_v,q) {
  let v = Math.abs(_v);
  let m = ((_v<0)? -1 : 1);
  let vq = Math.floor(v/q);
  if (v==0) { return 0; }
  if ((vq%2)==1) { return m*vq; }
  return m*(vq+1);
}

// floor version
//
//function _divq_floor(v,q) {
// let u = [];
//  for (let i=0; i<v.length; i++) {
//    u.push( Math.floor(v[i] / q) );
//  }
//  return u;
//}

// round to zero version
//
function _divq(v,q) {
  let u = [];
  for (let i=0; i<v.length; i++) {
    let _v = Math.abs(v[i]);
    let m = ((v[i]<0) ? -1 : 1);
    u.push( m*Math.floor(_v / q) );
  }
  return u;
}

function _div2(v) { return _divq(v,2); }

function _sgn(v) {
  if (v>0) { return  1; }
  if (v<0) { return -1; }
  return 0;
}

function _neg(v) {
  if (v.length==2) {
    return [-v[0], -v[1]];
  }
  return [ -v[0], -v[1], -v[2] ];
}

function _add(u,v) {
  if ((u.length == 2) || (v.length == 2)) {
    return [ u[0]+v[0], u[1]+v[1] ];
  }
  return [ u[0]+v[0], u[1]+v[1], u[2]+v[2] ];
}

function _mul(c,v) {
  if (v.length == 2) {
    return [ c*v[0], c*v[1] ];
  }
  return [ c*v[0], c*v[1], c*v[2] ];
}

function _dot(u,v) {
  if ((u.length == 2) || (v.length == 2)) {
    return  (u[0]*v[0]) + (u[1]*v[1]);
  }
  return (u[0]*v[0]) + (u[1]*v[1]) + (u[2]*v[2]);
}

function _abs(v) {
  let s = 0;
  for (let i=0; i<v.length; i++) {
    s += Math.abs(v[i]);
  }
  return s;
}

function _delta(v) {
  let u = [];
  for (let i=0; i<v.length; i++) {
    u.push( _sgn(v[i]) );
  }
  return u;
}

function _print(v) {
  if (v.length == 2)  { console.log(v[0], v[1]); }
  else                { console.log(v[0], v[1], v[2]); }
}

function _clone(v) {
  let u = [];
  for (let i=0; i<v.length; i++) {
    u.push(v[i]);
  }
  return u;
}

// Test to see if q is within bounds of volume
// whose corner is at p and volume defined by a,b,g
//
// If volume coordinate is positive, the corresponding p
// coordinate is taken to be lower bound.
// Otherwise, if the volume coordinate is negative,
// corresponding p coordinate is taken to be the
// upper bound.
//
// Works in 2 and 3 dimensions.
// Assumes q,p,a,b,g are all simple arrays (of length 2 or 3)
//
// q - query point
// p - corner point
// a - width like dimension
// b - height like dimension
// g - depth like dimension
//
function _inBounds(q, p, a, b, g) {
  let _a = [0,0,0],
      _b = [0,0,0],
      _g = [0,0,0];

  //let default_g = [0,0,1];
  //let _g = ((typeof g === "undefined") ? default_g : g);
  let _p = [0,0,0],
      _q = [0,0,0];

  _a[0] = a[0];
  _a[1] = a[1];
  _a[2] = ((a.length > 2) ? a[2] : 0);

  _b[0] = b[0];
  _b[1] = b[1];
  _b[2] = ((b.length > 2) ? b[2] : 0);

  _q[0] = q[0];
  _q[1] = q[1];
  _q[2] = ((q.length > 2) ? q[2] : 0);

  _p[0] = p[0];
  _p[1] = p[1];
  _p[2] = ((p.length > 2) ? p[2] : 0);

  if (typeof g === "undefined") {
    if      ((_a[0] == 0) && (_b[0] == 0)) { _g[0] = 1; }
    else if ((_a[1] == 0) && (_b[1] == 0)) { _g[1] = 1; }
    else if ((_a[2] == 0) && (_b[2] == 0)) { _g[2] = 1; }
  }
  else { _g = g; }


  let _d = [
    _a[0] + _b[0] + _g[0],
    _a[1] + _b[1] + _g[1],
    _a[2] + _b[2] + _g[2]
  ];


  for (let xyz=0; xyz<3; xyz++) {
    if ( _d[xyz] < 0 ) {
      if ((q[xyz] >  p[xyz]) ||
          (q[xyz] <= (p[xyz] + _d[xyz]))) { return false; }
    }
    else {
      if ((q[xyz] <  p[xyz]) ||
          (q[xyz] >= (p[xyz] + _d[xyz]))) { return false; }
    }
  }

  return true;

}

//-------------------
//         ___     __
//   ___ _|_  |___/ /
//  / _ `/ __// _  / 
//  \_, /____/\_,_/  
// /___/             
//-------------------


// "generalized" 2d gilbert curve
//
// alpha - width-like axis
// beta - height-like axis
//
// Enumerate points for the 2d Gilbert curve
// in alpha and beta axis.
// alpha/beta can be in 3d and should work properly.
//
// first prototype/reference implementation, single function, print only
//
function g2d_p(p, alpha, beta) {
  let a = _abs(alpha);
  let b = _abs(beta);

  let alpha2 = _div2(alpha);
  let beta2  = _div2(beta);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);

  let d_alpha = _delta(alpha);
  let d_beta  = _delta(beta);

  if (b==1) {

    let u = _clone(p);
    for (let i=0; i<a; i++) {
      _print(u);
      u = _add(u, d_alpha);
    }
    return;
  }

  if (a==1) {
    let u = _clone(p);
    for (let i=0; i<b; i++) {
      _print(u);
      u = _add(u, d_beta);
    }
    return;
  }

  if ( (2*a) > (3*b) ) {
    if ((a2%2) && (a>2)) { alpha2 = _add(alpha2, d_alpha); }


    g2d_p(p, alpha2, beta);
    g2d_p( _add(p, alpha2),
           _add(alpha, _neg(alpha2)),
           beta);

    return;
  }


  if ((b2%2) && (b>2)) { beta2 = _add(beta2, d_beta); }

  g2d_p( p,
         beta2,
         alpha2);
  g2d_p( _add(p, beta2),
         alpha,
         _add(beta, _neg(beta2)) );
  g2d_p( _add(p,
              _add( _add(alpha, _neg(d_alpha) ),
                    _add(beta2, _neg( d_beta) ) ) ),
         _neg(beta2),
         _add(alpha2, _neg(alpha)) );

}

//--------------------------------------------------------------
//   ______ ____           __  ___  ___       _____             
//  / ___(_) / /  ___ ____/ /_|_  |/ _ \  ___/ /_  |_ ____ _____
// / (_ / / / _ \/ -_) __/ __/ __// // / / _  / __/\ \ / // /_ /
// \___/_/_/_.__/\__/_/  \__/____/____/  \_,_/____/_\_\\_, //__/
//                                                    /___/     
//--------------------------------------------------------------

// "generalized" 2d gilbert curve
//
// alpha - width-like axis
// beta - height-like axis
//
// Enumerate points for the 2d Gilbert curve
// in alpha and beta axis.
// alpha/beta can be in 3d and should work properly.
//
// recursive, async
//
function Gilbert2D_d2xyz(dst_idx, cur_idx, p, alpha, beta) {
  let a = _abs(alpha);
  let b = _abs(beta);

  _dprint("#Gilbert2D_d2xyz: dst_idx:", dst_idx, "cur_idx:", cur_idx, "alpha:", alpha, "beta:", beta, "(", a, b, ")");

  let d_alpha = _delta(alpha);
  let d_beta  = _delta(beta);

  if (b==1) {

    _dprint("#Gilbert2D_d2xyz.Lb: dst_idx:", dst_idx, "cur_idx:", cur_idx, "alpha:", alpha, "beta:", beta, "(", a, b, ")");

    let d_idx = dst_idx - cur_idx;
    return _add(p, _mul(d_idx, d_alpha));
  }


  if (a==1) {

    _dprint("#Gilbert2D_d2xyz.La: dst_idx:", dst_idx, "cur_idx:", cur_idx, "alpha:", alpha, "beta:", beta, "(", a, b, ")");

    let d_idx = dst_idx - cur_idx;
    return _add(p, _mul(d_idx, d_beta));
  }

  let alpha2 = _div2(alpha);
  let beta2  = _div2(beta);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);

  _dprint("#  Gilbert2D_d2xyz: dst_idx:", dst_idx, "cur_idx:", cur_idx, "(alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

  if ( (2*a) > (3*b) ) {
    if ((a2%2) && (a>2)) {
      alpha2 = _add(alpha2, d_alpha);
      a2 = _abs(alpha2);
    }


    _dprint("#Gilbert2D_d2xyz.S0: dst_idx:", dst_idx, "cur_idx:", cur_idx, "(alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

    let nxt_idx = cur_idx + (a2*b);
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return Gilbert2D_d2xyz( dst_idx, cur_idx,
                              p,
                              alpha2,
                              beta );
    }
    cur_idx = nxt_idx;

    _dprint("#Gilbert2D_d2xyz.S1: dst_idx:", dst_idx, "cur_idx:", cur_idx, "(alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

    return Gilbert2D_d2xyz( dst_idx, cur_idx,
                            _add(p, alpha2),
                            _add(alpha, _neg(alpha2)),
                            beta );

  }


  if ((b2%2) && (b>2)) {
    beta2 = _add(beta2, d_beta);
    b2 = _abs(beta2);
  }


  _dprint("#Gilbert2D_d2xyz.A: dst_idx:", dst_idx, "cur_idx:", cur_idx, "(alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

  let nxt_idx = cur_idx + (b2*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert2D_d2xyz( dst_idx, cur_idx,
                            p,
                            beta2,
                            alpha2 );
  }
  cur_idx = nxt_idx;

  _dprint("#Gilbert2D_d2xyz.B: dst_idx:", dst_idx, "cur_idx:", cur_idx, "(alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

  nxt_idx = cur_idx + (a*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert2D_d2xyz( dst_idx, cur_idx,
                            _add(p, beta2),
                            alpha,
                            _add(beta, _neg(beta2)) );
  }
  cur_idx = nxt_idx;

  _dprint("#Gilbert2D_d2xyz.C: dst_idx:", dst_idx, "cur_idx:", cur_idx, "(alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

  return Gilbert2D_d2xyz( dst_idx, cur_idx,
                          _add(p,
                          _add( _add(alpha, _neg(d_alpha) ),
                                _add(beta2, _neg( d_beta) ) ) ),
                          _neg(beta2),
                          _add(alpha2, _neg(alpha)) );
}

//-----------------------------------------------------------------
//   ______ ____           __  ___  ___                  ___     __
//  / ___(_) / /  ___ ____/ /_|_  |/ _ \  __ ____ _____ |_  |___/ /
// / (_ / / / _ \/ -_) __/ __/ __// // /  \ \ / // /_ // __// _  / 
// \___/_/_/_.__/\__/_/  \__/____/____/  /_\_\\_, //__/____/\_,_/  
//                                           /___/                 
//-----------------------------------------------------------------

// "generalized" 2d gilbert curve
//
// alpha - width-like axis
// beta - height-like axis
//
// Enumerate points for the 2d Gilbert curve
// in alpha and beta axis.
// alpha/beta can be in 3d and should work properly.
//
// recursive, synchronous
//
function Gilbert2D_xyz2d(cur_idx, q, p, alpha, beta) {
  let a = _abs(alpha);
  let b = _abs(beta);

  let u = _clone(p);

  _dprint("#Gilbert2D_xyz2d: cur_idx:", cur_idx, "q:", q, "p:", p, "(alpha:", alpha, "beta:", beta, ") (", a, b, ")");

  let d_alpha = _delta(alpha);
  let d_beta  = _delta(beta);

  if ( b == 1 ) {
    return cur_idx + _dot( d_alpha, _add(q, _neg(u)) );
  }

  if ( a == 1 ) {
    return cur_idx + _dot( d_beta, _add(q, _neg(u)) );
  }

  let alpha2 = _div2(alpha);
  let beta2  = _div2(beta);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);

  _dprint("#  Gilbert2D_xyz2d: cur_idx:", cur_idx, "(alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

  if ( (2*a) > (3*b) ) {
    if ((a2%2) && (a>2)) {
      alpha2 = _add(alpha2, d_alpha);
      a2 = _abs(alpha2);
    }

    _dprint("#  Gilbert2D_xyz2d.W0");

    if (_inBounds(q, u, alpha2, beta)) {
      return Gilbert2D_xyz2d( cur_idx, q,
                              u,
                              alpha2,
                              beta );
    }
    cur_idx += (a2*b);
    u = _add(u, alpha2);

    _dprint("#  Gilbert2D_xyz2d.W1");

    return Gilbert2D_xyz2d( cur_idx, q,
                            u,
                            _add(alpha, _neg(alpha2)),
                            beta );

  }


  if ((b2%2) && (b>2)) {
    beta2 = _add(beta2, d_beta);
    b2 = _abs(beta2);
  }

  _dprint("#  Gilbert2D_xyz2d.A: q:", q, "u:", u, "beta2:", beta2, "alpha2:", alpha2);

  if (_inBounds(q, u, beta2, alpha2)) {
    return Gilbert2D_xyz2d( cur_idx, q,
                            u,
                            beta2,
                            alpha2 );
  }
  cur_idx += (b2*a2);

  _dprint("#  Gilbert2D_xyz2d.B");

  u = _add(p, beta2);
  if (_inBounds(q, u, alpha, _add(beta, _neg(beta2)))) {
    return Gilbert2D_xyz2d( cur_idx, q,
                            u,
                            alpha,
                            _add(beta, _neg(beta2)) );
  }
  cur_idx += (a*(b-b2));

  _dprint("#  Gilbert2D_xyz2d.C");

  u = _add(p,
           _add( _add(alpha, _neg(d_alpha) ),
                 _add(beta2, _neg( d_beta) ) ) );
  return Gilbert2D_xyz2d( cur_idx, q,
                          u,
                          _neg(beta2),
                          _add(alpha2, _neg(alpha)) );

}


//---------------------------------------------------------------
//   ______ ____           __  ___  ___  ___                    
//  / ___(_) / /  ___ ____/ /_|_  |/ _ \/ _ | ___ __ _____  ____
// / (_ / / / _ \/ -_) __/ __/ __// // / __ |(_-</ // / _ \/ __/
// \___/_/_/_.__/\__/_/  \__/____/____/_/ |_/___/\_, /_//_/\__/ 
//                                              /___/           
//---------------------------------------------------------------

// "generalized" 2d gilbert curve
//
// alpha - width-like axis
// beta - height-like axis
//
// Enumerate points for the 2d Gilbert curve
// in alpha and beta axis.
// alpha/beta can be in 3d and should work properly.
//
// recursive, async
//
function *Gilbert2DAsync(p, alpha, beta) {
  let a = _abs(alpha);
  let b = _abs(beta);

  _dprint("#Gilbert2DAsync", alpha, beta, "(", a, b, ")");

  let d_alpha = _delta(alpha);
  let d_beta  = _delta(beta);

  if (b==1) {
    let u = _clone(p);
    for (let i=0; i<a; i++) {
      yield u;
      u = _add(u, d_alpha);
    }
    return;
  }


  if (a==1) {
    let u = _clone(p);
    for (let i=0; i<b; i++) {
      yield u;
      u = _add(u, d_beta);
    }
    return;
  }

  let alpha2 = _div2(alpha);
  let beta2  = _div2(beta);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);

  _dprint("#  Gilbert2DAsync: (alpha2,beta2):", alpha2, beta2, "(", a2, b2, ")");

  if ( (2*a) > (3*b) ) {
    if ((a2%2) && (a>2)) { alpha2 = _add(alpha2, d_alpha); }

    _dprint("#  Gilbert2DAsync: 2a>3b:A");

    yield* Gilbert2DAsync( p, alpha2, beta );

    _dprint("#  Gilbert2DAsync: 2a>3b:B");

    yield* Gilbert2DAsync( _add(p, alpha2),
                           _add(alpha, _neg(alpha2)),
                           beta );

    return;
  }


  if ((b2%2) && (b>2)) { beta2 = _add(beta2, d_beta); }

  _dprint("#  Gilbert2DAsync: A (alpha:", alpha, "beta:", beta, "alpha2:", alpha2, "beta2:", beta2, ")", "(|a2|:", a2, "|b2|:", b2, ")");

  yield* Gilbert2DAsync( p,
                         beta2,
                         alpha2 );

  _dprint("#  Gilbert2DAsync: B (alpha:", alpha, "beta:", beta, "alpha2:", alpha2, "beta2:", beta2, ")", "(|a2|:", a2, "|b2|:", b2, ")");

  yield* Gilbert2DAsync( _add(p, beta2),
                         alpha,
                         _add(beta, _neg(beta2)) );

  _dprint("#  Gilbert2DAsync: C (alpha:", alpha, "beta:", beta, "alpha2:", alpha2, "beta2:", beta2, ")", "(|a2|:", a2, "|b2|:", b2, ")");

  yield* Gilbert2DAsync( _add(p,
                         _add( _add(alpha, _neg(d_alpha) ),
                               _add(beta2, _neg( d_beta) ) ) ),
                         _neg(beta2),
                         _add(alpha2, _neg(alpha)) );

  return;
}

//----------------------------------------
//                             ____    __
//   ___ __ _____  ____  ___ _|_  /___/ /
//  (_-</ // / _ \/ __/ / _ `//_ </ _  / 
// /___/\_, /_//_/\__/  \_, /____/\_,_/  
//     /___/           /___/             
//----------------------------------------

function Hilbert2x2x2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma) {
  let d_alpha = _delta(alpha);
  let d_beta  = _delta(beta);
  let d_gamma = _delta(gamma);

  let d_idx = dst_idx - cur_idx;

  _dprint("#H2x2x2: d_idx:", d_idx, "dst_idx:", dst_idx, "cur_idx:", cur_idx);

  switch (d_idx) {
    case 0: return [ p[0], p[1], p[2] ]; break;
    case 1: return _add(p, d_beta); break;
    case 2: return _add(p, _add(d_beta, d_gamma)); break;
    case 3: return _add(p, d_gamma); break;
    case 4: return _add(p, _add(d_alpha, d_gamma)); break;
    case 5: return _add(p, _add(d_alpha, _add(d_beta, d_gamma))); break;
    case 6: return _add(p, _add(d_alpha, d_beta)); break;
    case 7: return _add(p, d_alpha); break;
    default: return [-1,-1,-1]; break;
  }

  return [-1,-1,-1];
}

function Hilbert2x2x2_xyz2d(idx, q, p, alpha, beta, gamma) {
  let d_alpha = _delta(alpha);
  let d_beta  = _delta(beta);
  let d_gamma = _delta(gamma);

  let lu = [ 0, 7,  1, 6,  3, 4,  2, 5 ];

  let dxyz = [
    q[0] - p[0],
    q[1] - p[1],
    q[2] - p[2]
  ];

  let m_qp = _add(q, _neg(p));

  dxyz = [
    _dot(d_alpha, m_qp),
    _dot(d_beta, m_qp),
    _dot(d_gamma, m_qp)
  ];

  let p_idx = (4*dxyz[2]) + (2*dxyz[1]) + dxyz[0];

  _dprint("#Hilbert2x2x2_xyz2d: idx:", idx, "q:", q, "p:", p, "alpha:", alpha, "beta:", beta, "gamma:", gamma, "dxyz:", dxyz);

  if ((p_idx < 0) ||
      (p_idx > 7)) {
    return -1;
  }

  return idx + lu[p_idx];
}

//------------------------------------------------
//    ____    __  _______       _____             
//   |_  /___/ / / __/ _ \  ___/ /_  |_ ____ _____
//  _/_ </ _  / _\ \/ // / / _  / __/\ \ / // /_ /
// /____/\_,_/ /___/\___/  \_,_/____/_\_\\_, //__/
//                                      /___/     
//------------------------------------------------

function Gilbert3DS0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let d_alpha = _delta(alpha);

  let a   = _abs(alpha);
  let a2  = _abs(alpha2);

  let b   = _abs(beta);
  let g   = _abs(gamma);

  _dprint("#S0_d2xyz (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2)==1)) {
    alpha2 = _add(alpha2, d_alpha);
    a2 = _abs(alpha2);
  }

  let nxt_idx = cur_idx + (a2*b*g);

  _dprint("#S0.A_d2xyz {dst:",dst_idx, "cur:", cur_idx, "nxt:", nxt_idx, "}");

  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            p,
                            alpha2, beta, gamma );
  }
  cur_idx = nxt_idx;

  _dprint("#S0.B_d2xyz");

  return Gilbert3D_d2xyz( dst_idx, cur_idx,
                          _add(p, alpha2),
                          _add(alpha, _neg(alpha2)), beta, gamma );
}

//---------------------------------------------------
//    ____    __  _______                  ___     __
//   |_  /___/ / / __/ _ \  __ ____ _____ |_  |___/ /
//  _/_ </ _  / _\ \/ // /  \ \ / // /_ // __// _  / 
// /____/\_,_/ /___/\___/  /_\_\\_, //__/____/\_,_/  
//                             /___/                 
//---------------------------------------------------

function Gilbert3DS0_xyz2d( cur_idx, q, p, alpha, beta, gamma ) {

  let alpha2  = _div2(alpha);
  let d_alpha = _delta(alpha);

  let a   = _abs(alpha);
  let a2  = _abs(alpha2);

  let b   = _abs(beta);
  let g   = _abs(gamma);

  _dprint("#S0_xyz2d (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2)==1)) {
    alpha2 = _add(alpha2, d_alpha);
    a2 = _abs(alpha2);
  }

  _dprint("#S0.A_xyz2d");

  let u = _clone(p);
  if (_inBounds( q, p, alpha2, beta, gamma )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            alpha2, beta, gamma );
  }
  cur_idx += (a2*b*g);

  _dprint("#S0.B_xyz2d");

  u = _add(p, alpha2);
  return Gilbert3D_xyz2d( cur_idx, q,
                          u,
                          _add(alpha, _neg(alpha2)), beta, gamma );
}

//----------------------------------------------
//    ____    __  _______     _____             
//   |_  /___/ / / __<  / ___/ /_  |_ ____ _____
//  _/_ </ _  / _\ \ / / / _  / __/\ \ / // /_ /
// /____/\_,_/ /___//_/  \_,_/____/_\_\\_, //__/
//                                    /___/     
//----------------------------------------------

function Gilbert3DS1_d2xyz( dst_idx, cur_idx, p, alpha, beta, gamma ) {
  let alpha2 = _div2(alpha);
  let gamma3 = _divq(gamma, 3);

  let d_alpha = _delta(alpha);
  let d_gamma = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let g3 = _abs(gamma3);

  _dprint("#S1_d2xyz (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) {
    alpha2 = _add(alpha2, d_alpha);
    a2 = _abs(alpha2);
  }

  if ((g > 2) && ((g3 % 2) == 1)) {
    gamma3 = _add(gamma3, d_gamma);
    g3 = _abs(gamma3);
  }

  _dprint("#S1.A_d2xyz");

  let nxt_idx = cur_idx + (g3*a2*b);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            p,
                            gamma3, alpha2, beta );
  }
  cur_idx = nxt_idx;

  _dprint("#S1.B_d2xyz");

  nxt_idx = cur_idx + (a*b*(g-g3));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add(p, gamma3),
                            alpha, beta, _add(gamma, _neg(gamma3)) );
  }
  cur_idx = nxt_idx;

  _dprint("#S1.C_d2xyz");

  return Gilbert3D_d2xyz( dst_idx, cur_idx,
                          _add(p, _add( _add(alpha, _neg(d_alpha)), _add(gamma3, _neg(d_gamma)) ) ),
                          _neg(gamma3), _neg(_add(alpha, _neg(alpha2))), beta );
}

//-------------------------------------------------
//    ____    __  _______                ___     __
//   |_  /___/ / / __<  / __ ____ _____ |_  |___/ /
//  _/_ </ _  / _\ \ / /  \ \ / // /_ // __// _  / 
// /____/\_,_/ /___//_/  /_\_\\_, //__/____/\_,_/  
//                           /___/                 
//-------------------------------------------------

function Gilbert3DS1_xyz2d( cur_idx, q, p, alpha, beta, gamma ) {
  let alpha2 = _div2(alpha);
  let gamma3 = _divq(gamma, 3);

  let d_alpha = _delta(alpha);
  let d_gamma = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let g3 = _abs(gamma3);

  _dprint("#S1_xyz2d (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) {
    alpha2 = _add(alpha2, d_alpha);
    a2 = _abs(alpha2);
  }

  if ((g > 2) && ((g3 % 2) == 1)) {
    gamma3 = _add(gamma3, d_gamma);
    g3 = _abs(gamma3);
  }

  _dprint("#S1.A_xyz2d");

  let u = _clone(p);
  if (_inBounds(q, u, gamma3, alpha2, beta )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            gamma3, alpha2, beta );
  }
  cur_idx += (g3*a2*b);

  _dprint("#S1.B_xyz2d");

  u =_add(p, gamma3);
  if (_inBounds(q, u, alpha, beta, _add(gamma, _neg(gamma3)) )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            alpha, beta, _add(gamma, _neg(gamma3)) );

  }
  cur_idx += (a*b*(g-g3));

  _dprint("#S1.C_xyz2d");

  u = _add(p, _add( _add(alpha, _neg(d_alpha)), _add(gamma3, _neg(d_gamma)) ) );
  return Gilbert3D_xyz2d( cur_idx, q,
                          u,
                          _neg(gamma3), _neg(_add(alpha, _neg(alpha2))), beta );
}



//------------------------------------------------
//    ____    __  _______       _____             
//   |_  /___/ / / __/_  |  ___/ /_  |_ ____ _____
//  _/_ </ _  / _\ \/ __/  / _  / __/\ \ / // /_ /
// /____/\_,_/ /___/____/  \_,_/____/_\_\\_, //__/
//                                      /___/     
//------------------------------------------------

function Gilbert3DS2_d2xyz( dst_idx, cur_idx, p, alpha, beta, gamma ) {
  let alpha2 = _div2(alpha);
  let beta3 = _divq(beta, 3);

  let d_alpha = _delta(alpha);
  let d_beta = _delta(beta);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b3 = _abs(beta3);

  _dprint("#S2_d2xyz (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) {
    alpha2 = _add(alpha2, d_alpha);
    a2 = _abs(alpha2);
  }

  if ((b > 2) && ((b3 % 2) == 1)) {
    beta3 = _add(beta3, d_beta);
    b3 = _abs(beta3);
  }

  _dprint("#S2.A_d2xyz");

  let nxt_idx = cur_idx + (b3*g*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            p,
                            beta3, gamma, alpha2 );
  }
  cur_idx = nxt_idx;

  _dprint("#S2.B_d2xyz");

  nxt_idx = cur_idx + (a*(b-b3)*g);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add(p, beta3),
                            alpha, _add(beta, _neg(beta3)), gamma );

  }
  cur_idx = nxt_idx;

  _dprint("#S2.C_d2xyz");

  return Gilbert3D_d2xyz( dst_idx, cur_idx,
                          _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta3, _neg(d_beta)) ) ),
                          _neg(beta3), gamma, _neg(_add(alpha, _neg(alpha2))) );

}


//----------------------------------------------------
//    ____    __  _______                  ___     __
//   |_  /___/ / / __/_  |  __ ____ _____ |_  |___/ /
//  _/_ </ _  / _\ \/ __/   \ \ / // /_ // __// _  / 
// /____/\_,_/ /___/____/  /_\_\\_, //__/____/\_,_/  
//                             /___/                 
//----------------------------------------------------

function Gilbert3DS2_xyz2d(cur_idx, q, p, alpha, beta, gamma) {
  let alpha2 = _div2(alpha);
  let beta3 = _divq(beta, 3);

  let d_alpha = _delta(alpha);
  let d_beta = _delta(beta);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b3 = _abs(beta3);

  _dprint("#S2_xyz2d cur_idx:", cur_idx, "q:", q, "p:", p, "(a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) {
    alpha2 = _add(alpha2, d_alpha);
    a2 = _abs(alpha2);
  }

  if ((b > 2) && ((b3 % 2) == 1)) {
    beta3 = _add(beta3, d_beta);
    b3 = _abs(beta3);
  }

  _dprint("#S2.A_xyz2d");

  let u = _clone(p);
  if (_inBounds(q, u, beta3, gamma, alpha2 )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            beta3, gamma, alpha2 );
  }
  cur_idx += (b3*g*a2);

  _dprint("#S2.B_xyz2d");

  u = _add(p, beta3);
  if (_inBounds(q, u, alpha, _add(beta, _neg(beta3)), gamma )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            alpha, _add(beta, _neg(beta3)), gamma );
  }
  cur_idx += (a*(b-b3)*g);

  _dprint("#S2.C_xyz2d");

  u = _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta3, _neg(d_beta)) ) );
  return Gilbert3D_xyz2d( cur_idx, q,
                          u,
                          _neg(beta3), gamma, _neg(_add(alpha, _neg(alpha2))) );

}


//-------------------------------------------------
//    ____    __     _____       _____             
//   |_  /___/ / __ / / _ \  ___/ /_  |_ ____ _____
//  _/_ </ _  / / // / // / / _  / __/\ \ / // /_ /
// /____/\_,_/  \___/\___/  \_,_/____/_\_\\_, //__/
//                                       /___/     
//-------------------------------------------------

function Gilbert3DJ0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma) {
  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J0_d2xyz (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) { alpha2 = _add(alpha2, d_alpha); a2 = _abs(alpha2); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta);   b2 = _abs(beta2); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); g2 = _abs(gamma2); }

  _dprint("#J0.A_d2xyz");

  let nxt_idx = cur_idx + (b2*g2*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            p,
                            beta2, gamma2, alpha2 );
  }
  cur_idx = nxt_idx;

  _dprint("#J0.B_d2xyz");

  nxt_idx = cur_idx + (g*a2*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add(p, beta2),
                            gamma, alpha2, _add(beta, _neg(beta2)) );
  }
  cur_idx = nxt_idx;

  _dprint("#J0.C_d2xyz");

  nxt_idx = cur_idx + (a*b2*(g-g2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add( p, _add( _add(beta2, _neg(d_beta)), _add(gamma, _neg(d_gamma)) ) ),
                            alpha, _neg(beta2), _neg( _add(gamma, _neg(gamma2)) ) );
  }
  cur_idx = nxt_idx;

  _dprint("#J0.D_d2xyz");

  nxt_idx = cur_idx + (g*(a-a2)*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta2, _add(gamma, _neg(d_gamma))) ) ),
                            _neg(gamma), _neg( _add(alpha, _neg(alpha2)) ), _add(beta, _neg(beta2)) );
  }
  cur_idx = nxt_idx;

  _dprint("#J0.E_d2xyz");

  return Gilbert3D_d2xyz( dst_idx, cur_idx,
                          _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta2, _neg(d_beta)) ) ),
                          _neg(beta2), gamma2, _neg( _add(alpha, _neg(alpha2)) ) );
}

//----------------------------------------------------
//    ____    __     _____                  ___     __
//   |_  /___/ / __ / / _ \  __ ____ _____ |_  |___/ /
//  _/_ </ _  / / // / // /  \ \ / // /_ // __// _  / 
// /____/\_,_/  \___/\___/  /_\_\\_, //__/____/\_,_/  
//                              /___/                 
//----------------------------------------------------

function Gilbert3DJ0_xyz2d(cur_idx, q, p, alpha, beta, gamma) {
  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J0_xyz2d (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) { alpha2 = _add(alpha2, d_alpha); a2 = _abs(alpha2); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta);   b2 = _abs(beta2); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); g2 = _abs(gamma2); }

  _dprint("#J0.A_xyz2d");

  let u = _clone(p);
  if (_inBounds(q, u, beta2, gamma2, alpha2)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            beta2, gamma2, alpha2 );
  }
  cur_idx += (b2*g2*a2);

  _dprint("#J0.B_xyz2d");

  u = _add(p, beta2);
  if (_inBounds(q, u, gamma, alpha2, _add(beta, _neg(beta2)) )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            gamma, alpha2, _add(beta, _neg(beta2)) );
  }
  cur_idx += (g*a2*(b-b2));

  _dprint("#J0.C_xyz2d");

  u = _add( p, _add( _add(beta2, _neg(d_beta)), _add(gamma, _neg(d_gamma)) ) );
  if (_inBounds(q, u, alpha, _neg(beta2), _neg( _add(gamma, _neg(gamma2)) ) )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            alpha, _neg(beta2), _neg( _add(gamma, _neg(gamma2)) ) );
  }
  cur_idx += (a*b2*(g-g2));

  _dprint("#J0.D_xyz2d");

  u = _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta2, _add(gamma, _neg(d_gamma))) ) );
  if (_inBounds(q, u, _neg(gamma), _neg( _add(alpha, _neg(alpha2)) ), _add(beta, _neg(beta2)) )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            _neg(gamma), _neg( _add(alpha, _neg(alpha2)) ), _add(beta, _neg(beta2)) );
  }
  cur_idx += (g*(a-a2)*(b-b2));

  _dprint("#J0.E_xyz2d");

  u = _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta2, _neg(d_beta)) ) );
  return Gilbert3D_xyz2d( cur_idx, q,
                          u,
                          _neg(beta2), gamma2, _neg( _add(alpha, _neg(alpha2)) ) );
}

//-----------------------------------------------
//    ____    __     _____     _____             
//   |_  /___/ / __ / <  / ___/ /_  |_ ____ _____
//  _/_ </ _  / / // // / / _  / __/\ \ / // /_ /
// /____/\_,_/  \___//_/  \_,_/____/_\_\\_, //__/
//                                     /___/     
//-----------------------------------------------

function Gilbert3DJ1_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J1_d2xyz (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 0)) { alpha2 = _add(alpha2, d_alpha); a2 = _abs(alpha2); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta);   b2 = _abs(beta2); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); g2 = _abs(gamma2); }

  _dprint("#J1.A_d2xyz");

  let nxt_idx = cur_idx + (g2*a2*b2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            p,
                            gamma2, alpha2, beta2 );
  }
  cur_idx = nxt_idx;

  _dprint("#J1.B_d2xyz");

  nxt_idx = cur_idx + (b*(g-g2)*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add( p, gamma2 ),
                            beta, _add(gamma, _neg(gamma2)), alpha2 );
  }
  cur_idx = nxt_idx;

  _dprint("#J1.C_d2xyz");

  nxt_idx = cur_idx + (a*(b-b2)*g2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add( p, _add( _add(gamma2, _neg(d_gamma)), _add(beta, _neg(d_beta)) ) ),
                            alpha, _neg(_add(beta, _neg(beta2))), _neg(gamma2) );
  }
  cur_idx = nxt_idx;

  _dprint("#J1.D_d2xyz");

  nxt_idx = cur_idx + (b*(g-g2)*(a-a2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add( p , _add( _add(alpha, _neg(d_alpha)), _add( _add(beta, _neg(d_beta)), gamma2 ) ) ),
                            _neg(beta), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))) );
  }
  cur_idx = nxt_idx;

  _dprint("#J1.E_d2xyz");

  return Gilbert3D_d2xyz( dst_idx, cur_idx,
                          _add( p, _add( _add(alpha, _neg(d_alpha)), _add(gamma2, _neg(d_gamma)) ) ),
                          _neg(gamma2), _neg(_add(alpha, _neg(alpha2))), beta2 );

}

//--------------------------------------------------
//    ____    __     _____                ___     __
//   |_  /___/ / __ / <  / __ ____ _____ |_  |___/ /
//  _/_ </ _  / / // // /  \ \ / // /_ // __// _  / 
// /____/\_,_/  \___//_/  /_\_\\_, //__/____/\_,_/  
//                            /___/                 
//--------------------------------------------------

function Gilbert3DJ1_xyz2d(cur_idx, q, p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J1_xyz2d (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 0)) { alpha2 = _add(alpha2, d_alpha); a2 = _abs(alpha2); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta);   b2 = _abs(beta2); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); g2 = _abs(gamma2); }

  _dprint("#J1.A_xyz2d");

  let u = _clone(p);
  if (_inBounds(q, u, gamma2, alpha2, beta2)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            p,
                            gamma2, alpha2, beta2 );
  }
  cur_idx += (g2*a2*b2);

  _dprint("#J1.B_xyz2d");

  u = _add(p, gamma2);
  if (_inBounds(q, u, beta, _add(gamma, _neg(gamma2)), alpha2)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            beta, _add(gamma, _neg(gamma2)), alpha2 );
  }
  cur_idx += b*(g-g2)*a2;

  _dprint("#J1.C_xyz2d");

  u = _add( p, _add( _add(gamma2, _neg(d_gamma)), _add(beta, _neg(d_beta)) ) );
  if (_inBounds(q, u, alpha, _neg(_add(beta, _neg(beta2))), _neg(gamma2) )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            alpha, _neg(_add(beta, _neg(beta2))), _neg(gamma2) );
  }
  cur_idx += a*(b-b2)*g2;

  _dprint("#J1.D_xyz2d");

  u = _add( p , _add( _add(alpha, _neg(d_alpha)), _add( _add(beta, _neg(d_beta)), gamma2 ) ) );
  if (_inBounds(q, u, _neg(beta), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))) )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            _neg(beta), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))) );
  }
  cur_idx += (b*(g-g2)*(a-a2));

  _dprint("#J1.E_xyz2d");

  u = _add( p, _add( _add(alpha, _neg(d_alpha)), _add(gamma2, _neg(d_gamma)) ) );
  return Gilbert3D_xyz2d( cur_idx, q,
                          u,
                          _neg(gamma2), _neg(_add(alpha, _neg(alpha2))), beta2 );

}

//-----------------------------------------------
//    ____    __     _____       _____             
//   |_  /___/ / __ / /_  |  ___/ /_  |_ ____ _____
//  _/_ </ _  / / // / __/  / _  / __/\ \ / // /_ /
// /____/\_,_/  \___/____/  \_,_/____/_\_\\_, //__/
//                                       /___/     
//-----------------------------------------------

function Gilbert3DJ2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J2_d2xyz (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 0)) { alpha2 = _add(alpha2, d_alpha); a2 = _abs(alpha2); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta);   b2 = _abs(beta2); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); g2 = _abs(gamma2); }

  _dprint("#J2.A_d2xyz");

  let nxt_idx = cur_idx + (b2*g*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            p,
                            beta2, gamma, alpha2 );
  }
  cur_idx = nxt_idx;

  _dprint("#J2.B_d2xyz");

  nxt_idx = cur_idx + (g2*a*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add(p, beta2),
                            gamma2, alpha, _add(beta, _neg(beta2)) );
  }
  cur_idx = nxt_idx;

  _dprint("#J2.C_d2xyz");

  nxt_idx = cur_idx + (a*(b-b2)*(g-g2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add(p, _add(beta2, gamma2)),
                            alpha, _add(beta, _neg(beta2)), _add(gamma, _neg(gamma2)) );
  }
  cur_idx = nxt_idx;

  _dprint("#J2.D_d2xyz");

  nxt_idx = cur_idx + (b2*(g-g2)*(a-a2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            _add( p, _add( _add( alpha, _neg(d_alpha) ), _add( _add(beta2, _neg(d_beta)), gamma2 ) ) ),
                            _neg(beta2), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))) );
  }
  cur_idx = nxt_idx;

  _dprint("#J2.E_d2xyz");

  return Gilbert3D_d2xyz( dst_idx, cur_idx,
                          _add( p, _add( _add(alpha, _neg(d_alpha)), _add( gamma2, _neg(d_gamma) ) ) ),
                          _neg(gamma2), _neg(_add(alpha, _neg(alpha2))), beta2);

}

//----------------------------------------------------
//    ____    __     _____                  ___     __
//   |_  /___/ / __ / /_  |  __ ____ _____ |_  |___/ /
//  _/_ </ _  / / // / __/   \ \ / // /_ // __// _  / 
// /____/\_,_/  \___/____/  /_\_\\_, //__/____/\_,_/  
//                              /___/                 
//----------------------------------------------------


function Gilbert3DJ2_xyz2d(cur_idx, q, p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);


  _dprint("#J2_xyz2d (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 0)) { alpha2 = _add(alpha2, d_alpha); a2 = _abs(alpha2); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta);   b2 = _abs(beta2); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); g2 = _abs(gamma2); }

  _dprint("#J2.A_xyz2d");

  let u = _clone(p);
  if (_inBounds(q, p, beta2, gamma, alpha2)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            beta2, gamma, alpha2 );
  }
  cur_idx += (b2*g*a2);

  _dprint("#J2.B_xyz2d");

  u = _add(p, beta2);
  if (_inBounds(q, u, gamma2, alpha, _add(beta, _neg(beta2)))) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            gamma2, alpha, _add(beta, _neg(beta2)) );
  }
  cur_idx += (g2*a*(b-b2));

  _dprint("#J2.C_xyz2d");

  u = _add(p, _add(beta2, gamma2));
  if (_inBounds(q, u, alpha, _add(beta, _neg(beta2)), _add(gamma, _neg(gamma2)))) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            alpha, _add(beta, _neg(beta2)), _add(gamma, _neg(gamma2)) );
  }
  cur_idx += (a*(b-b2)*(g-g2));

  _dprint("#J2.D_xyz2d");

  u = _add( p, _add( _add( alpha, _neg(d_alpha) ), _add( _add(beta2, _neg(d_beta)), gamma2 ) ) );
  if (_inBounds(q, u, _neg(beta2), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))))) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            u,
                            _neg(beta2), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))) );
  }
  cur_idx += (b2*(g-g2)*(a-a2));

  _dprint("#J2.E_xyz2d");

  u = _add( p, _add( _add(alpha, _neg(d_alpha)), _add( gamma2, _neg(d_gamma) ) ) );
  return Gilbert3D_xyz2d( cur_idx, q,
                          u,
                          _neg(gamma2), _neg(_add(alpha, _neg(alpha2))), beta2);

}

//--------------------------------------------------------------
//   ______ ____           __  ____ ___       _____             
//  / ___(_) / /  ___ ____/ /_|_  // _ \  ___/ /_  |_ ____ _____
// / (_ / / / _ \/ -_) __/ __//_ </ // / / _  / __/\ \ / // /_ /
// \___/_/_/_.__/\__/_/  \__/____/____/  \_,_/____/_\_\\_, //__/
//                                                    /___/     
//--------------------------------------------------------------

// Gilbert3d d2xyz
//
function Gilbert3D_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma) {
  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a0 = (a % 2);
  let b0 = (b % 2);
  let g0 = (g % 2);

  _dprint("#Gilbert3D_d2xyz: dst_idx:", dst_idx, "cur_idx:", cur_idx, "p:", p, "a:", alpha, "b:", beta, "g:", gamma);

  // base cases
  //
  if ((a == 2) &&
      (b == 2) &&
      (g == 2)) {

    _dprint("#H2x2x2_d2xyz:");

    return Hilbert2x2x2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if (a == 1) { return Gilbert2D_d2xyz(dst_idx, cur_idx, p, beta, gamma); }
  if (b == 1) { return Gilbert2D_d2xyz(dst_idx, cur_idx, p, alpha, gamma); }
  if (g == 1) { return Gilbert2D_d2xyz(dst_idx, cur_idx, p, alpha, beta); }

  // eccentric cases
  //
  if (((3*a) > (5*b)) &&
      ((3*a) > (5*g))) {

    _dprint("#S0_d2xyz:");

    return Gilbert3DS0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if (((2*b) > (3*g)) ||
      ((2*b) > (3*a))) {

    _dprint("#S2_d2xyz:");

    return Gilbert3DS2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if ((2*g) > (3*b)) {

    _dprint("#S1_d2xyz:");

    return Gilbert3DS1_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  // bulk recursion
  //
  if (g0 == 0) {

    _dprint("#J0_d2xyz:");

    return Gilbert3DJ0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if ((a0 == 0) || (b0 == 0)) {

    _dprint("#J1_d2xyz:");

    return Gilbert3DJ1_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  _dprint("#J2_d2xyz:");

  // a0 == b0 == g0 == 1
  //
  return Gilbert3DJ2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma);
}

//-----------------------------------------------------------------
//   ______ ____           __  ____ ___                  ___     __
//  / ___(_) / /  ___ ____/ /_|_  // _ \  __ ____ _____ |_  |___/ /
// / (_ / / / _ \/ -_) __/ __//_ </ // /  \ \ / // /_ // __// _  / 
// \___/_/_/_.__/\__/_/  \__/____/____/  /_\_\\_, //__/____/\_,_/  
//                                           /___/                 
//-----------------------------------------------------------------

// Gilbert3d xyz2d
//
function Gilbert3D_xyz2d(cur_idx, q, p, alpha, beta, gamma) {
  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a0 = (a % 2);
  let b0 = (b % 2);
  let g0 = (g % 2);

  _dprint("#Gilbert3D_xyz2d: cur_idx:", cur_idx, "p:", p, "a:", alpha, "b:", beta, "g:", gamma);

  // base cases
  //
  if ((a == 2) &&
      (b == 2) &&
      (g == 2)) {

    _dprint("#H2x2x2_xyz2d:");

    return Hilbert2x2x2_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if (a == 1) { return Gilbert2D_xyz2d(cur_idx, q, p, beta, gamma); }
  if (b == 1) { return Gilbert2D_xyz2d(cur_idx, q, p, alpha, gamma); }
  if (g == 1) { return Gilbert2D_xyz2d(cur_idx, q, p, alpha, beta); }

  // eccentric cases
  //
  if (((3*a) > (5*b)) &&
      ((3*a) > (5*g))) {

    _dprint("#S0_xyz2d:");

    return Gilbert3DS0_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if (((2*b) > (3*g)) ||
      ((2*b) > (3*a))) {

    _dprint("#S2_xyz2d:");

    return Gilbert3DS2_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if ((2*g) > (3*b)) {

    _dprint("#S1_xyz2d:");

    return Gilbert3DS1_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  // bulk recursion
  //
  if (g0 == 0) {

    _dprint("#J0_xyz2d:");

    return Gilbert3DJ0_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if ((a0 == 0) || (b0 == 0)) {

    _dprint("#J1_xyz2d:");

    return Gilbert3DJ1_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  _dprint("#J2_xyz2d:");

  // a0 == b0 == g0 == 1
  //
  return Gilbert3DJ2_xyz2d(cur_idx, q, p, alpha, beta, gamma);
}


//--------------------------------------------
//                                 ____    __
//  ___ ____ __ _____  ____  ___ _|_  /___/ /
// / _ `(_-</ // / _ \/ __/ / _ `//_ </ _  / 
// \_,_/___/\_, /_//_/\__/  \_, /____/\_,_/  
//         /___/           /___/             
//--------------------------------------------


function *Hilbert2x2x2Async(p, alpha, beta, gamma) {
  let xyz = _clone(p);

  let d_alpha = _delta(alpha);
  let d_beta  = _delta(beta);
  let d_gamma = _delta(gamma);

  yield [ xyz[0], xyz[1], xyz[2] ];

  xyz = _add(xyz, d_beta);
  yield [ xyz[0], xyz[1], xyz[2] ];

  xyz = _add(xyz, d_gamma);
  yield [ xyz[0], xyz[1], xyz[2] ];

  xyz = _add(xyz, _neg(d_beta));
  yield [ xyz[0], xyz[1], xyz[2] ];

  xyz = _add(xyz, d_alpha);
  yield [ xyz[0], xyz[1], xyz[2] ];

  xyz = _add(xyz, d_beta);
  yield [ xyz[0], xyz[1], xyz[2] ];

  xyz = _add(xyz, _neg(d_gamma));
  yield [ xyz[0], xyz[1], xyz[2] ];

  xyz = _add(xyz, _neg(d_beta));
  yield [ xyz[0], xyz[1], xyz[2] ];
}

function *Gilbert3DS0Async(p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let d_alpha = _delta(alpha);

  let a   = _abs(alpha);
  let a2  = _abs(alpha2);

  _dprint("#S0 (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2)==1)) {
    alpha2 = _add(alpha2, d_alpha);
  }

  _dprint("#S0.A");

  yield *Gilbert3DAsync( p,
                         alpha2, beta, gamma );

  _dprint("#S0.B");

  yield *Gilbert3DAsync( _add(p, alpha2),
                         _add(alpha, _neg(alpha2)), beta, gamma );
}

function *Gilbert3DS1Async(p, alpha, beta, gamma) {
  let alpha2 = _div2(alpha);
  let gamma3 = _divq(gamma, 3);

  let d_alpha = _delta(alpha);
  let d_gamma = _delta(gamma);

  let a = _abs(alpha);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let g3 = _abs(gamma3);

  _dprint("#S1 (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) {
    alpha2 = _add(alpha2, d_alpha);
  }

  if ((g > 2) && ((g3 % 2) == 1)) {
    gamma3 = _add(gamma3, d_gamma);
  }

  _dprint("#S1.A");

  yield *Gilbert3DAsync( p,
                         gamma3, alpha2, beta );

  _dprint("#S1.B");

  yield *Gilbert3DAsync( _add(p, gamma3),
                         alpha, beta, _add(gamma, _neg(gamma3)) );

  _dprint("#S1.C");

  yield *Gilbert3DAsync( _add(p, _add( _add(alpha, _neg(d_alpha)), _add(gamma3, _neg(d_gamma)) ) ),
                         _neg(gamma3), _neg(_add(alpha, _neg(alpha2))), beta );
}

function *Gilbert3DS2Async(p, alpha, beta, gamma) {
  let alpha2 = _div2(alpha);
  let beta3 = _divq(beta, 3);

  let d_alpha = _delta(alpha);
  let d_beta = _delta(beta);

  let a = _abs(alpha);
  let b = _abs(beta);

  let a2 = _abs(alpha2);
  let b3 = _abs(beta3);

  _dprint("#S2 (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) {
    alpha2 = _add(alpha2, d_alpha);
  }

  if ((b > 2) && ((b3 % 2) == 1)) {
    beta3 = _add(beta3, d_beta);
  }

  _dprint("#S2.A");

  yield *Gilbert3DAsync( p,
                         beta3, gamma, alpha2 );

  _dprint("#S2.B");

  yield *Gilbert3DAsync( _add(p, beta3),
                         alpha, _add(beta, _neg(beta3)), gamma );

  _dprint("#S2.C");

  yield *Gilbert3DAsync( _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta3, _neg(d_beta)) ) ),
                         _neg(beta3), gamma, _neg(_add(alpha, _neg(alpha2))) );

}

function *Gilbert3DJ0Async(p, alpha, beta, gamma) {
  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J0 (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 1)) { alpha2 = _add(alpha2, d_alpha); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); }

  _dprint("#J0.A");

  yield *Gilbert3DAsync( p,
                         beta2, gamma2, alpha2 );

  _dprint("#J0.B");

  yield *Gilbert3DAsync( _add(p, beta2),
                         gamma, alpha2, _add(beta, _neg(beta2)) );

  _dprint("#J0.C");

  yield *Gilbert3DAsync( _add( p, _add( _add(beta2, _neg(d_beta)), _add(gamma, _neg(d_gamma)) ) ),
                         alpha, _neg(beta2), _neg( _add(gamma, _neg(gamma2)) ) );

  _dprint("#J0.D");

  yield *Gilbert3DAsync( _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta2, _add(gamma, _neg(d_gamma))) ) ),
                         _neg(gamma), _neg( _add(alpha, _neg(alpha2)) ), _add(beta, _neg(beta2)) );

  _dprint("#J0.E");

  yield *Gilbert3DAsync( _add( p, _add( _add(alpha, _neg(d_alpha)), _add(beta2, _neg(d_beta)) ) ),
                         _neg(beta2), gamma2, _neg( _add(alpha, _neg(alpha2)) ) );
}

function *Gilbert3DJ1Async(p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J1 (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 0)) { alpha2 = _add(alpha2, d_alpha); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); }

  _dprint("#J1.A");

  yield *Gilbert3DAsync(p,
                        gamma2, alpha2, beta2 );


  _dprint("#J1.B");

  yield *Gilbert3DAsync(_add( p, gamma2 ),
                        beta, _add(gamma, _neg(gamma2)), alpha2 );

  _dprint("#J1.C");

  yield *Gilbert3DAsync(_add( p, _add( _add(gamma2, _neg(d_gamma)), _add(beta, _neg(d_beta)) ) ),
                        alpha, _neg(_add(beta, _neg(beta2))), _neg(gamma2) );

  _dprint("#J1.D");

  yield *Gilbert3DAsync(_add( p , _add( _add(alpha, _neg(d_alpha)), _add( _add(beta, _neg(d_beta)), gamma2 ) ) ),
                        _neg(beta), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))) );

  _dprint("#J1.E");

  yield *Gilbert3DAsync(_add( p, _add( _add(alpha, _neg(d_alpha)), _add(gamma2, _neg(d_gamma)) ) ),
                        _neg(gamma2), _neg(_add(alpha, _neg(alpha2))), beta2 );

}

function *Gilbert3DJ2Async(p, alpha, beta, gamma) {

  let alpha2  = _div2(alpha);
  let beta2   = _div2(beta);
  let gamma2  = _div2(gamma);

  let d_alpha  = _delta(alpha);
  let d_beta   = _delta(beta);
  let d_gamma  = _delta(gamma);

  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a2 = _abs(alpha2);
  let b2 = _abs(beta2);
  let g2 = _abs(gamma2);

  _dprint("#J2 (a:", _abs(alpha), "b:", _abs(beta), "g:", _abs(gamma), ")");

  if ((a > 2) && ((a2 % 2) == 0)) { alpha2 = _add(alpha2, d_alpha); }
  if ((b > 2) && ((b2 % 2) == 1)) { beta2  = _add(beta2, d_beta); }
  if ((g > 2) && ((g2 % 2) == 1)) { gamma2 = _add(gamma2, d_gamma); }

  _dprint("#J2.A");

  yield *Gilbert3DAsync( p,
                         beta2, gamma, alpha2 );

  _dprint("#J2.B");

  yield *Gilbert3DAsync( _add(p, beta2),
                         gamma2, alpha, _add(beta, _neg(beta2)) );

  _dprint("#J2.C");

  yield *Gilbert3DAsync( _add(p, _add(beta2, gamma2)),
                         alpha, _add(beta, _neg(beta2)), _add(gamma, _neg(gamma2)) );

  _dprint("#J2.D");

  yield *Gilbert3DAsync( _add( p, _add( _add( alpha, _neg(d_alpha) ), _add( _add(beta2, _neg(d_beta)), gamma2 ) ) ),
                         _neg(beta2), _add(gamma, _neg(gamma2)), _neg(_add(alpha, _neg(alpha2))) );

  _dprint("#J2.E");

  yield *Gilbert3DAsync( _add( p, _add( _add(alpha, _neg(d_alpha)), _add( gamma2, _neg(d_gamma) ) ) ),
                         _neg(gamma2), _neg(_add(alpha, _neg(alpha2))), beta2);

}

// Gilbert3dAsync
//
function *Gilbert3DAsync(p, alpha, beta, gamma) {
  let a = _abs(alpha);
  let b = _abs(beta);
  let g = _abs(gamma);

  let a0 = (a % 2);
  let b0 = (b % 2);
  let g0 = (g % 2);

  _dprint("#Gilbert3DAsync: p:", p, "a:", alpha, "b:", beta, "g:", gamma);

  // base cases
  //
  if ((a == 2) &&
      (b == 2) &&
      (g == 2)) {

    _dprint("#H2x2x2:");

    yield *Hilbert2x2x2Async(p, alpha, beta, gamma);
    return;
  }

  if (a == 1) { yield *Gilbert2DAsync(p, beta, gamma); return; }
  if (b == 1) { yield *Gilbert2DAsync(p, alpha, gamma); return; }
  if (g == 1) { yield *Gilbert2DAsync(p, alpha, beta); return; }

  // eccentric cases
  //
  if (((3*a) > (5*b)) &&
      ((3*a) > (5*g))) {

    _dprint("#S0:");

    yield *Gilbert3DS0Async(p, alpha, beta, gamma);
    return;
  }

  if (((2*b) > (3*g)) ||
      ((2*b) > (3*a))) {

    _dprint("#S2:");

    yield *Gilbert3DS2Async(p, alpha, beta, gamma);
    return;
  }

  if ((2*g) > (3*b)) {

    _dprint("#S1:");

    yield *Gilbert3DS1Async(p, alpha, beta, gamma);
    return;
  }

  // bulk recursion
  //
  if (g0 == 0) {

    _dprint("#J0:");

    yield *Gilbert3DJ0Async(p, alpha, beta, gamma);
    return;
  }

  if ((a0 == 0) || (b0 == 0)) {

    _dprint("#J1:");

    yield *Gilbert3DJ1Async(p, alpha, beta, gamma);
    return;
  }

  _dprint("#J2:");

  // a0 == b0 == g0 == 1
  //
  yield *Gilbert3DJ2Async(p, alpha, beta, gamma);
}


//--------------------------------------------

//---
//---
//---

//---
//---
//---

// Gilbert3D++
//
function Gilbert3D(w, h, d) {
  let p = [0,0,0],
      alpha = [w,0,0],
      beta = [0,h,0],
      gamma = [0,0,d];

  let pnt = [];

  let g3xyz = Gilbert3DAsync(p, alpha, beta, gamma);
  for (let hv = g3xyz.next() ; !hv.done ; hv = g3xyz.next()) {
    let v = hv.value;
    pnt.push( [v[0], v[1], v[2]] );
  }

  return pnt;
}

function Gilbert2D(w,h) {
  let p = [0,0,0],
      alpha = [w,0,0],
      beta = [0,h,0];

  let pnt = [];

  let g2xy = Gilbert2DAsync(p, alpha, beta);
  for (let hv = g2xy.next() ; !hv.done ; hv = g2xy.next()) {
    let v = hv.value;
    pnt.push( [v[0], v[1]] );
  }

  return pnt;
}

var OP_LIST = [
  "xy", "xyz", "xyzp",
  "xy2d", "d2xy",
  "xyz2d", "d2xyz",

  "xyz2da", "d2xyza",
  "xy2da", "d2xya"
];


function _show_help(msg) {
  msg = ((typeof msg !== "undefined") ? msg : "");
  if (msg.length > 0) { console.log(msg, "\n"); }

  let op_list = OP_LIST;

  console.log("");
  console.log("usage:");
  console.log("");
  console.log("  node ./gilbert3dpp.js [OP] [W] [H] [D]");
  console.log("");
  console.log("  OP   one of '" + op_list.join("', '") + "'.");
  console.log("       For adapative methods (xy2da,d2xya,xyz2da,d2xyza), the suffix '.0', '.1' or '.2'");
  console.log("       can be added to indicate which adapted method to use (0:harmony,1:hamiltonian,2:axis).");
  console.log("  W    width");
  console.log("  H    height");
  console.log("  D    depth");
  console.log("");
  console.log("example:");
  console.log("");
  console.log("  node ./gilbert3dpp.js d2xya.1 10 12");
  console.log("");
}

function _main(argv) {
  let op = "xyz";
  let w = 0;
  let h = 0;
  let d = 1;

  let arg_idx = 1;

  let op_list = OP_LIST;

  if (argv.length <= 1) {
    _show_help();
  }
  else {

    if (argv.length > 1) { op = argv[1]; }
    if (argv.length > 2) { w = parseInt(argv[2]); }
    if (argv.length > 3) { h = parseInt(argv[3]); }
    if (argv.length > 4) { d = parseInt(argv[4]); }

    let adapt_method = GILBERT_ADAPT_METHOD.HARMONY;

    let op_tok = op.split(".");
    if (op_tok.length > 1) {
      op = op_tok[0];
      adapt_method = parseInt(op_tok[1]);
    }

    let op_found = false
    for (let i=0; i<op_list.length; i++) {
      if (op == op_list[i]) { op_found = true; break; }
    }
    if (!op_found) {
      _show_help();
      return;
    }

    if (isNaN(w) || (w==0) ||
        isNaN(h) || (h==0) ||
        isNaN(d) || (d==0)) {
      _show_help();
      return;
    }

    if (op == "xy") {
      let p = Gilbert2D(w,h);
      for (let i=0; i<p.length; i++) {
        console.log(p[i][0], p[i][1]);
      }
    }
    else if (op == "xyz") {
      let p = Gilbert3D(w,h,d);
      for (let i=0; i<p.length; i++) {
        console.log(p[i][0], p[i][1], p[i][2]);
      }
    }

    else if (op == "xy2d") {

      for (let y=0; y<h; y++) {
        for (let x=0; x<w; x++) {
          let idx = Gilbert2D_xyz2d(0, [x,y,0], [0,0,0], [w,0,0], [0,h,0]);
          console.log(idx, x,y);
        }
      }

    }

    else if (op == "d2xy") {

      for (let idx=0; idx<(w*h); idx++) {
        let xyz = Gilbert2D_d2xyz(idx, 0, [0,0,0], [w,0,0], [0,h,0]);
        console.log(xyz[0], xyz[1]);
      }

    }

    else if (op == "xyz2d") {

      for (let z=0; z<d; z++) {
        for (let y=0; y<h; y++) {
          for (let x=0; x<w; x++) {
            let idx = Gilbert3D_xyz2d(0, [x,y,z], [0,0,0], [w,0,0], [0,h,0], [0,0,d]);
            console.log(idx, x,y,z);
          }
        }
      }

    }

    else if (op == "d2xyz") {

      for (let idx=0; idx<(w*h*d); idx++) {
        let xyz = Gilbert3D_d2xyz(idx, 0, [0,0,0], [w,0,0], [0,h,0], [0,0,d]);
        console.log(xyz[0], xyz[1], xyz[2]);
      }

    }


    // "dynamic" functions that adjust the endpoints to create a Hamiltonian path
    //

    else if (op == "d2xyza") {

      for (let idx=0; idx<(w*h*d); idx++) {
        let xyz = Gilbert3DAdapt_d2xyz(idx, w,h,d, adapt_method);
        console.log(xyz[0], xyz[1], xyz[2]);
      }

    }

    else if (op == "xyz2da") {

      for (let z=0; z<d; z++) {
        for (let y=0; y<h; y++) {
          for (let x=0; x<w; x++) {
            let idx = Gilbert3DAdapt_xyz2d([x,y,z], w,h,d, adapt_method);
            console.log(idx, x,y,z);
          }
        }
      }

    }

    else if (op == "d2xya") {

      for (let idx=0; idx<(w*h*d); idx++) {
        let xyz = Gilbert2DAdapt_d2xy(idx, w,h, adapt_method);
        console.log(xyz[0], xyz[1]);
      }

    }

    else if (op == "xy2da") {

      for (let y=0; y<h; y++) {
        for (let x=0; x<w; x++) {
          let idx = Gilbert2DAdapt_xy2d([x,y], w,h, adapt_method);
          console.log(idx, x,y);
        }
      }

    }

  }
}

if (typeof module !== "undefined") {
  module.exports["Gilbert2D"] = Gilbert2D;
  module.exports["Gilbert3D"] = Gilbert3D;

  module.exports["Gilbert2DAsync"] = Gilbert2DAsync;
  module.exports["Gilbert3DAsync"] = Gilbert3DAsync;

  module.exports["Gilbert2D_d2xyz"] = Gilbert2D_d2xyz;
  module.exports["Gilbert2D_xyz2d"] = Gilbert2D_xyz2d;

  module.exports["Gilbert3D_d2xyz"] = Gilbert3D_d2xyz;
  module.exports["Gilbert3D_xyz2d"] = Gilbert3D_xyz2d;

  module.exports["Gilbert2DAdapt_d2xyz"] = Gilbert2DAdapt_d2xy;
  module.exports["Gilbert2DAdapt_xyz2d"] = Gilbert2DAdapt_xy2d;

  module.exports["Gilbert3DAdapt_d2xyz"] = Gilbert3DAdapt_d2xyz;
  module.exports["Gilbert3DAdapt_xyz2d"] = Gilbert3DAdapt_xyz2d;

  module.exports["ADAPT_METHOD"] = GILBERT_ADAPT_METHOD;
}

//---
// see https://stackoverflow.com/questions/4981891/node-js-equivalent-of-pythons-if-name-main
//
if ((typeof require !== "undefined")  &&
    (require.main === module)) {
  _main(process.argv.slice(1));
}
//---

