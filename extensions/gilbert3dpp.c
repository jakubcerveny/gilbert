/*

 To the extent possible under law, the person who associated CC0 with
 this project has waived all copyright and related or neighboring rights
 to this project.

 You should have received a copy of the CC0 legalcode along with this
 work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

*/

#include <stdio.h>
#include <stdlib.h>

enum GILBERT3DPP_ADAPT_METHOD {
  HARMONY = 0,
  HAMILTONIAN,
  AXIS
};

/*****************************************************************
 *    __       __               ___              __  _
 *   / /  ___ / /__  ___ ____  / _/_ _____  ____/ /_(_)__  ___  ___
 *  / _ \/ -_) / _ \/ -_) __/ / _/ // / _ \/ __/ __/ / _ \/ _ \(_-<
 * /_//_/\__/_/ .__/\__/_/   /_/ \_,_/_//_/\__/\__/_/\___/_//_/___/
 *           /_/
 ******************************************************************/

int d2e(int _v) {
  int v, m, v2;

  v = abs(_v);
  m = ((_v < 0) ? -1 : 1);
  v2 = v/2;

  if (v == 0) { return 0; }
  if ((v2%2) == 0) { return m*v2; }
  return m*(v2+1);
}

int d2u(int _v) {
  int v, m, v2;

  v = abs(_v);
  m = ((_v < 0) ? -1 : 1);
  v2 = (v/2);

  if (v==0) { return 0; }
  if ((v2%2) == 1) { return m*v2; }
  return m*(v2+1);
}

int dqe(int _v,int q) {
  int v, m, vq;

  v = abs(_v);
  m = ((_v < 0) ? -1 : 1);
  vq = (v/q);

  if ((vq%2) == 0) { return m*vq; }
  return m*(vq+1);
}

int dqu(int _v, int q) {
  int v, m, vq;

  v = abs(_v);
  m = ((_v < 0) ? -1 : 1);
  vq = (v/q);

  if (v==0) { return 0; }
  if ((vq%2)==1) { return m*vq; }
  return m*(vq+1);
}

int v_divq(int *u, int *v, int q) {
  int i, _v, m;

  for (i=0; i<3; i++) {
    _v = abs(v[i]);
    m = ( (v[i]<0) ? -1 : 1 );
    u[i] = m*(_v / q);
  }

  return 0;
}

int v_div2(int *u, int *v) { return v_divq(u,v,2); }

int sgn(int v) {
  if (v < 0) { return -1; }
  if (v > 0) { return  1; }
  return 0;
}

int v_neg(int *u, int *v) {
  u[0] = -v[0];
  u[1] = -v[1];
  u[2] = -v[2];
  return 0;
}

int v_add(int *u, int *v, int *w) {
  u[0] = v[0] + w[0];
  u[1] = v[1] + w[1];
  u[2] = v[2] + w[2];
}

int v_sub(int *u, int *v, int *w) {
  u[0] = v[0] - w[0];
  u[1] = v[1] - w[1];
  u[2] = v[2] - w[2];
}

int v_mul(int *u, int c, int *v) {
  u[0] = c*v[0];
  u[1] = c*v[1];
  u[2] = c*v[2];
}

int dot_v(int *v, int *w) {
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

int abs_sum_v(int *v) {
  return abs(v[0]) + abs(v[1]) + abs(v[2]);
}

int v_delta(int *u, int *v) {
  u[0] = sgn(v[0]);
  u[1] = sgn(v[1]);
  u[2] = sgn(v[2]);
  return 0;
}

int inBounds(int *q, int *p, int *alpha, int *beta, int *_gamma) {
  int xyz, d;
  int gamma[3] = {0};

  if (!_gamma) {
    if      ( (alpha[0] == 0) && (beta[0] == 0) ) { gamma[0] = 1; }
    else if ( (alpha[1] == 0) && (beta[1] == 0) ) { gamma[1] = 1; }
    else if ( (alpha[2] == 0) && (beta[2] == 0) ) { gamma[2] = 1; }
  }
  else {
    gamma[0] = _gamma[0];
    gamma[1] = _gamma[1];
    gamma[2] = _gamma[2];
  }


  for (xyz=0; xyz<3; xyz++) {
    d = alpha[xyz] + beta[xyz] + gamma[xyz];

    if (d < 0) {
      if ((q[xyz] > p[xyz]) ||
          (q[xyz] <= (p[xyz] + d))) { return 0; }
    }
    else {
      if ((q[xyz] < p[xyz]) ||
          (q[xyz] >= (p[xyz] + d))) { return 0; }
    }

  }

  return 1;
}

int inBounds2(int *q, int *p, int *alpha, int *beta) {
  return inBounds(q, p, alpha, beta, NULL);
}

int Hilbert2x2x2_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma) {
  int d_alpha[3], d_beta[3], d_gamma[3],
      tv[3];
  int d_idx;

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  d_idx = dst_idx - cur_idx;

  u[0] = -1;
  u[1] = -1;
  u[2] = -1;

  switch (d_idx) {
    case 0:
      u[0]=p[0];
      u[1]=p[1];
      u[2]=p[2];
      break;
    case 1:
      v_add(u, p, d_beta);
      break;
    case 2:
      v_add(tv, d_beta, d_gamma);
      v_add(u, p, tv);
      break;
    case 3:
      v_add(u, p, d_gamma);
      break;
    case 4:
      v_add(tv, d_alpha, d_gamma);
      v_add(u, p, tv);
      break;
    case 5:
      v_add(tv, d_beta, d_gamma);
      v_add(tv, tv, d_alpha);
      v_add(u, p, tv);
      break;
    case 6:
      v_add(tv, d_alpha, d_beta);
      v_add(u, p, tv);
      break;
    case 7:
      v_add(u, p, d_alpha);
      break;
    default:
      return -1;
      break;
  }

  return 0;
}

int Hilbert2x2x2_xyz2d(int idx, int *q, int *p, int *alpha, int *beta, int *gamma) {
  int d_alpha[3], d_beta[3], d_gamma[3],
      dxyz[3], m_qp[3];
  int lu[8] = {0, 7,  1, 6,  3, 4,  2, 5 };
  int p_idx=0;

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  v_sub(m_qp, q,p);

  dxyz[0] = dot_v(d_alpha, m_qp);
  dxyz[1] = dot_v(d_beta, m_qp);
  dxyz[2] = dot_v(d_gamma, m_qp);

  p_idx = (4*dxyz[2]) + (2*dxyz[1]) + dxyz[0];

  if ((p_idx < 0) ||
      (p_idx > 7)) {
    return -1;
  }

  return idx + lu[p_idx];
}


/*****************************************************************
 *    __       __               ___              __  _
 *   / /  ___ / /__  ___ ____  / _/_ _____  ____/ /_(_)__  ___  ___
 *  / _ \/ -_) / _ \/ -_) __/ / _/ // / _ \/ __/ __/ / _ \/ _ \(_-<
 * /_//_/\__/_/ .__/\__/_/   /_/ \_,_/_//_/\__/\__/_/\___/_//_/___/
 *           /_/
 ******************************************************************/


/**************************************************
 *   ______ ____           __  ___  ___    __    __ 
 *  / ___(_) / /  ___ ____/ /_|_  |/ _ \__/ /___/ /_
 * / (_ / / / _ \/ -_) __/ __/ __// // /_  __/_  __/
 * \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
 *                                                  
 **************************************************/

/************************
 *     _ ___             
 *  __| |_  )_ ___  _ ___
 * / _` |/ /\ \ / || |_ /
 * \__,_/___/_\_\\_, /__|
 *               |__/    
 *************************/

/*
 * Index to xyz position
 * Generalized so that this works on 2d (axis-aligned) planes in 3d
 *
 *  u       - output vector (where converted answer is stored) (int[3])
 *  dst_idx - query index
 *  cur_idx - current index (pass in 0 if calling from outside)
 *  p       - start point (int[3])
 *  alpha   - width-like vector (int[3])
 *  beta    - height-like vector (int[3])
 *
 * Return:
 *
 *   0 - success
 *  !0 - failure
 *
 * Example:
 *
 *   int u_result[3],
 *       p[3] = {0},
 *       alpha[3] = {20,0,0},
 *       beta[3] = {0,44,0;
 *   Gilbert2D_d2xyz(u_result, 50, 0, p, alpha, beta);
 *   printf("index:%i maps to (%i,%i) (in width %i, height %i)"\n,
 *           50, u_result[0], u_result[1], 20, 44);
 *
 */
int Gilbert2D_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta) {
  int a, b,
      a2, b2,
      d_idx = 0;
  int alpha2[3], beta2[3],
      d_alpha[3], d_beta[3];
  int tv[3];

  int t_alpha[3], t_beta[3], t_p[3];

  int nxt_idx;

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);

  if (b==1) {
    d_idx = dst_idx - cur_idx;
    v_mul(tv, d_idx, d_alpha);
    v_add(u, p, tv);
    return 0;
  }


  if (a==1) {
    d_idx = dst_idx - cur_idx;
    v_mul(tv, d_idx, d_beta);
    v_add(u, p, tv);
    return 0;
  }

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);

  if ( (2*a) > (3*b) ) {
    if ((a2%2) && (a>2)) {
      v_add(alpha2, alpha2, d_alpha);
      a2 = abs_sum_v(alpha2);
    }

    nxt_idx = cur_idx + (a2*b);
    if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
      return Gilbert2D_d2xyz( u, dst_idx, cur_idx,
                              p,
                              alpha2,
                              beta );
    }
    cur_idx = nxt_idx;


    v_add(t_p, p, alpha2);
    v_sub(t_alpha, alpha, alpha2);

    return Gilbert2D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            t_alpha,
                            beta );

  }


  if ((b2%2) && (b>2)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }


  nxt_idx = cur_idx + (b2*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert2D_d2xyz( u, dst_idx, cur_idx,
                            p,
                            beta2,
                            alpha2 );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (a*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    v_add(t_p, p, beta2);
    v_sub(t_beta, beta, beta2);
    return Gilbert2D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            alpha,
                            t_beta );
  }
  cur_idx = nxt_idx;

  v_add(t_p, p, alpha);
  v_sub(t_p, t_p, d_alpha);

  v_add(t_p, t_p, beta2);
  v_sub(t_p, t_p, d_beta);

  v_neg(t_alpha, beta2);

  v_sub(t_beta, alpha2, alpha);

  return Gilbert2D_d2xyz( u, dst_idx, cur_idx,
                          t_p,
                          t_alpha,
                          t_beta );
}


/************************
 *              ___    _
 * __ ___  _ __|_  )__| |
 * \ \ / || |_ // // _` |
 * /_\_\\_, /__/___\__,_|
 *      |__/
 ************************/

/*
 * xyz position to index
 * Generalized so that this works on 2d (axis-aligned) planes in 3d
 *
 *  cur_idx - current index (pass in 0 if calling from outside)
 *  q       - query point (int[3])
 *  p       - start point (int[3])
 *  alpha   - width-like vector (int[3])
 *  beta    - height-like vector (int[3])
 *
 * Return:
 *
 *   0 - success
 *  !0 - failure
 *
 * Example:
 *
 *   int idx,
 *       q[3] = {8,4,0},
 *       p[3] = {0},
 *       alpha[3] = {20,0,0},
 *       beta[3] = {0,44,0;
 *   idx = Gilbert2D_xyz2d(0, q, p, alpha, beta);
 *   printf("point (%i,%i) maps to index %i (in width %i, height %i)"\n,
 *           q[0], q[1], idx, 20, 44);
 *
 */
int Gilbert2D_xyz2d(int cur_idx, int *q, int *p, int *alpha, int *beta) {
  int a,b, a2,b2;
  int d_alpha[3], d_beta[3],
      t_alpha[3], t_beta[3],
      alpha2[3], beta2[3];

  int tv[3], u[3];

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);       
                            
  u[0] = p[0];
  u[1] = p[1];
  u[2] = p[2];

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
    
  if ( b == 1 ) {
    v_sub(tv, q, u);
    return cur_idx + dot_v( d_alpha, tv );
  }

  if ( a == 1 ) {
    v_sub(tv, q, u);
    return cur_idx + dot_v( d_beta, tv );
  }

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);

  if ( (2*a) > (3*b) ) {
    if ((a2%2) && (a>2)) {
      v_add(alpha2, alpha2, d_alpha);
      a2 = abs_sum_v(alpha2);
    }

    if (inBounds2(q, u, alpha2, beta)) {
      return Gilbert2D_xyz2d( cur_idx, q,
                              u,
                              alpha2,
                              beta );
    }
    cur_idx += (a2*b);
    v_add(u, u, alpha2);

    v_sub(t_alpha, alpha, alpha2);
    return Gilbert2D_xyz2d( cur_idx, q,
                            u,
                            t_alpha,
                            beta );

  }


  if ((b2%2) && (b>2)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }

  if (inBounds2(q, u, beta2, alpha2)) {
    return Gilbert2D_xyz2d( cur_idx, q,
                            u,
                            beta2,
                            alpha2 );
  }
  cur_idx += (b2*a2);

  v_add(u, p, beta2);
  v_sub(t_beta, beta, beta2);
  if (inBounds2(q, u, alpha, t_beta)) {
    return Gilbert2D_xyz2d( cur_idx, q,
                            u,
                            alpha,
                            t_beta );
  }
  cur_idx += (a*(b-b2));

  v_add(u, p, alpha);
  v_sub(u, u, d_alpha);

  v_add(u, u, beta2);
  v_sub(u, u, d_beta);

  v_neg(t_alpha, beta2);
  v_sub(t_beta, alpha2, alpha);
  return Gilbert2D_xyz2d( cur_idx, q,
                          u,
                          t_alpha,
                          t_beta );

}

/***************************************************
 *   ______ ____           __  ___  ___    __    __ 
 *  / ___(_) / /  ___ ____/ /_|_  |/ _ \__/ /___/ /_
 * / (_ / / / _ \/ -_) __/ __/ __// // /_  __/_  __/
 * \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
 *                                                  
 **************************************************/


/***************************************************
 *   ______ ____           __  ____ ___    __    __ 
 *  / ___(_) / /  ___ ____/ /_|_  // _ \__/ /___/ /_
 * / (_ / / / _ \/ -_) __/ __//_ </ // /_  __/_  __/
 * \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
 *                                                  
 ***************************************************/

int Gilbert3D_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma);
int Gilbert3D_xyz2d(int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma);

int Gilbert3DS0_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma) {
  int alpha2[3], d_alpha[3], t_alpha[3];
  int a, a2,
      b, g;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_delta(d_alpha, alpha);

  a = abs_sum_v(alpha);
  a2 = abs_sum_v(alpha2);

  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  if ((a > 2) && ((a2 % 2)==1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  nxt_idx = cur_idx + (a2*b*g);

  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            p,
                            alpha2, beta, gamma );
  }
  cur_idx = nxt_idx;

  v_add(t_p, p, alpha2);
  v_sub(t_alpha, alpha, alpha2);
  return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                          t_p,
                          t_alpha, beta, gamma );
}

//---------------------------------------------------
//    ____    __  _______                  ___     __
//   |_  /___/ / / __/ _ \  __ ____ _____ |_  |___/ /
//  _/_ </ _  / _\ \/ // /  \ \ / // /_ // __// _  / 
// /____/\_,_/ /___/\___/  /_\_\\_, //__/____/\_,_/  
//                             /___/                 
//---------------------------------------------------

int Gilbert3DS0_xyz2d(int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma) {
  int alpha2[3], d_alpha[3], t_alpha[3];
  int a, a2;
  int b, g;
  int t_p[3];

  v_div2(alpha2, alpha);
  v_delta(d_alpha, alpha);

  a   = abs_sum_v(alpha);
  a2  = abs_sum_v(alpha2);

  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  if ((a > 2) && ((a2 % 2)==1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  t_p[0] = p[0];
  t_p[1] = p[1];
  t_p[2] = p[2];

  if (inBounds( q, p, alpha2, beta, gamma )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            alpha2, beta, gamma );
  }
  cur_idx += (a2*b*g);

  v_add(t_p, p, alpha2);
  v_sub(t_alpha, alpha, alpha2);
  return Gilbert3D_xyz2d( cur_idx, q,
                          t_p,
                          t_alpha, beta, gamma );
}

//----------------------------------------------
//    ____    __  _______     _____             
//   |_  /___/ / / __<  / ___/ /_  |_ ____ _____
//  _/_ </ _  / _\ \ / / / _  / __/\ \ / // /_ /
// /____/\_,_/ /___//_/  \_,_/____/_\_\\_, //__/
//                                    /___/     
//----------------------------------------------

int Gilbert3DS1_d2xyz( int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma ) {

  int alpha2[3], d_alpha[3], t_alpha[3],
      gamma3[3], d_gamma[3], t_gamma[3],
      t_beta[3];
  int a, b, g,
      a2, g3;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_divq(gamma3, gamma, 3);

  v_delta(d_alpha, alpha);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  g3 = abs_sum_v(gamma3);

  if ((a > 2) && ((a2 % 2) == 1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((g > 2) && ((g3 % 2) == 1)) {
    v_add(gamma3, gamma3, d_gamma);
    g3 = abs_sum_v(gamma3);
  }

  nxt_idx = cur_idx + (g3*a2*b);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            p,
                            gamma3, alpha2, beta );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (a*b*(g-g3));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    v_add(t_p, p, gamma3);
    v_sub(t_gamma, gamma, gamma3);
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            alpha, beta, t_gamma );
  }
  cur_idx = nxt_idx;

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, gamma3);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, gamma3);

  v_sub(t_beta, alpha2, alpha);
  return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                          t_p,
                          t_alpha, t_beta, beta );
}

//-------------------------------------------------
//    ____    __  _______                ___     __
//   |_  /___/ / / __<  / __ ____ _____ |_  |___/ /
//  _/_ </ _  / _\ \ / /  \ \ / // /_ // __// _  / 
// /____/\_,_/ /___//_/  /_\_\\_, //__/____/\_,_/  
//                           /___/                 
//-------------------------------------------------

int Gilbert3DS1_xyz2d( int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma ) {
  int alpha2[3], d_alpha[3], t_alpha[3],
      gamma3[3], d_gamma[3], t_gamma[3],
      t_beta[3];
  int a, b, g,
      a2, g3;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_divq(gamma3, gamma, 3);

  v_delta(d_alpha, alpha);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  g3 = abs_sum_v(gamma3);

  if ((a > 2) && ((a2 % 2) == 1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((g > 2) && ((g3 % 2) == 1)) {
    v_add(gamma3, gamma3, d_gamma);
    g3 = abs_sum_v(gamma3);
  }

  t_p[0] = p[0];
  t_p[1] = p[1];
  t_p[2] = p[2];

  if (inBounds( q, t_p, gamma3, alpha2, beta )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            gamma3, alpha2, beta );
  }
  cur_idx += (g3*a2*b);

  v_add(t_p, p, gamma3);
  v_sub(t_gamma, gamma, gamma3);
  if (inBounds(q, t_p, alpha, beta, t_gamma )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            alpha, beta, t_gamma );
  }
  cur_idx += (a*b*(g-g3));


  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, gamma3);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, gamma3);

  v_sub(t_beta, alpha2, alpha);

  return Gilbert3D_xyz2d( cur_idx, q,
                          t_p,
                          t_alpha, t_beta, beta );
}



//------------------------------------------------
//    ____    __  _______       _____             
//   |_  /___/ / / __/_  |  ___/ /_  |_ ____ _____
//  _/_ </ _  / _\ \/ __/  / _  / __/\ \ / // /_ /
// /____/\_,_/ /___/____/  \_,_/____/_\_\\_, //__/
//                                      /___/     
//------------------------------------------------

int Gilbert3DS2_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma) {

  int alpha2[3], d_alpha[3], t_alpha[3],
      beta3[3], d_beta[3], t_beta[3],
      t_gamma[3];
  int a, b, g,
      a2, b3;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_divq(beta3, beta, 3);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b3 = abs_sum_v(beta3);

  if ((a > 2) && ((a2 % 2) == 1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b3 % 2) == 1)) {
    v_add(beta3, beta3, d_beta);
    b3 = abs_sum_v(beta3);
  }

  nxt_idx = cur_idx + (b3*g*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            p,
                            beta3, gamma, alpha2 );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (a*(b-b3)*g);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    v_add(t_p, p, beta3);

    v_sub(t_beta, beta, beta3);

    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            alpha, t_beta, gamma );

  }
  cur_idx = nxt_idx;

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, beta3);
  v_sub(t_p, t_p, d_beta);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, beta3);
  v_sub(t_gamma, alpha2, alpha);

  return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                          t_p,
                          t_alpha, gamma, t_gamma );

}


//----------------------------------------------------
//    ____    __  _______                  ___     __
//   |_  /___/ / / __/_  |  __ ____ _____ |_  |___/ /
//  _/_ </ _  / _\ \/ __/   \ \ / // /_ // __// _  / 
// /____/\_,_/ /___/____/  /_\_\\_, //__/____/\_,_/  
//                             /___/                 
//----------------------------------------------------

int Gilbert3DS2_xyz2d( int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma ) {

  int alpha2[3], d_alpha[3], t_alpha[3],
      beta3[3], d_beta[3], t_beta[3],
      t_gamma[3];
  int a, b, g,
      a2, b3;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_divq(beta3, beta, 3);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b3 = abs_sum_v(beta3);

  if ((a > 2) && ((a2 % 2) == 1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b3 % 2) == 1)) {
    v_add(beta3, beta3, d_beta);
    b3 = abs_sum_v(beta3);
  }

  t_p[0] = p[0];
  t_p[1] = p[1];
  t_p[2] = p[2];

  if (inBounds( q, t_p, beta3, gamma, alpha2 )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            beta3, gamma, alpha2 );
  }
  cur_idx += (b3*g*a2);

  v_add(t_p, p, beta3);
  v_sub(t_beta, beta, beta3);
  if (inBounds(q, t_p, alpha, t_beta, gamma )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            alpha, t_beta, gamma );
  }
  cur_idx += (a*(b-b3)*g);

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, beta3);
  v_sub(t_p, t_p, d_beta);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, beta3);

  v_sub(t_gamma, alpha2, alpha);

  return Gilbert3D_xyz2d( cur_idx, q,
                          t_p,
                          t_alpha, gamma, t_gamma );
}


//-------------------------------------------------
//    ____    __     _____       _____             
//   |_  /___/ / __ / / _ \  ___/ /_  |_ ____ _____
//  _/_ </ _  / / // / // / / _  / __/\ \ / // /_ /
// /____/\_,_/  \___/\___/  \_,_/____/_\_\\_, //__/
//                                       /___/     
//-------------------------------------------------

int Gilbert3DJ0_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma) {


  int alpha2[3], d_alpha[3], t_alpha[3],
      beta2[3], d_beta[3], t_beta[3],
      gamma2[3], d_gamma[3], t_gamma[3];

  int a, b, g,
      a2, b2, g2;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);
  v_div2(gamma2, gamma);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);
  g2 = abs_sum_v(gamma2);

  if ((a > 2) && ((a2 % 2) == 1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b2 % 2) == 1)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }

  if ((g > 2) && ((g2 % 2) == 1)) {
    v_add(gamma2, gamma2, d_gamma);
    g2 = abs_sum_v(gamma2);
  }

  nxt_idx = cur_idx + (b2*g2*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            p,
                            beta2, gamma2, alpha2 );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (g*a2*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    v_add(t_p, p, beta2);
    v_sub(t_gamma, beta, beta2);
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            gamma, alpha2, t_gamma );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (a*b2*(g-g2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    v_sub(t_p, beta2, d_beta);
    v_add(t_p, t_p, gamma);
    v_sub(t_p, t_p, d_gamma);
    v_add(t_p, t_p, p);

    v_neg(t_beta, beta2);

    v_sub(t_gamma, gamma2, gamma);

    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            alpha, t_beta, t_gamma );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (g*(a-a2)*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    v_sub(t_p, alpha, d_alpha);
    v_add(t_p, t_p, beta2);
    v_add(t_p, t_p, gamma);
    v_sub(t_p, t_p, d_gamma);
    v_add(t_p, t_p, p);

    v_neg(t_alpha, gamma);
    
    v_sub(t_beta, alpha2, alpha);

    v_sub(t_gamma, beta, beta2);

    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            t_alpha, t_beta, t_gamma );
  }
  cur_idx = nxt_idx;

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, beta2);
  v_sub(t_p, t_p, d_beta);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, beta2);

  v_sub(t_gamma, alpha2, alpha);

  return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                          t_p,
                          t_alpha, gamma2, t_gamma );
}

//----------------------------------------------------
//    ____    __     _____                  ___     __
//   |_  /___/ / __ / / _ \  __ ____ _____ |_  |___/ /
//  _/_ </ _  / / // / // /  \ \ / // /_ // __// _  / 
// /____/\_,_/  \___/\___/  /_\_\\_, //__/____/\_,_/  
//                              /___/                 
//----------------------------------------------------

int Gilbert3DJ0_xyz2d(int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma) {

  int alpha2[3],  d_alpha[3], t_alpha[3],
      beta2[3],   d_beta[3],  t_beta[3],
      gamma2[3],  d_gamma[3], t_gamma[3];

  int a,  b,  g,
      a2, b2, g2;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);
  v_div2(gamma2, gamma);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);
  g2 = abs_sum_v(gamma2);

  if ((a > 2) && ((a2 % 2) == 1)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b2 % 2) == 1)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }

  if ((g > 2) && ((g2 % 2) == 1)) {
    v_add(gamma2, gamma2, d_gamma);
    g2 = abs_sum_v(gamma2);
  }

  t_p[0] = p[0];
  t_p[1] = p[1];
  t_p[2] = p[2];

  if (inBounds( q, t_p, beta2, gamma2, alpha2 )) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            beta2, gamma2, alpha2 );
  }
  cur_idx += (b2*g2*a2);

  v_add(t_p, p, beta2);
  v_sub(t_gamma, beta, beta2);
  if (inBounds(q, t_p, gamma, alpha2, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            gamma, alpha2, t_gamma );
  }
  cur_idx += (g*a2*(b-b2));

  v_sub(t_p, beta2, d_beta);
  v_add(t_p, t_p, gamma);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_beta, beta2);

  v_sub(t_gamma, gamma2, gamma);

  if (inBounds(q, t_p, alpha, t_beta, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            alpha, t_beta, t_gamma );
  }
  cur_idx += (a*b2*(g-g2));

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, beta2);
  v_add(t_p, t_p, gamma);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, gamma);

  v_sub(t_beta, alpha2, alpha);

  v_sub(t_gamma, beta, beta2);

  if (inBounds(q, t_p, t_alpha, t_beta, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            t_alpha, t_beta, t_gamma);
  }
  cur_idx += (g*(a-a2)*(b-b2));

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, beta2);
  v_sub(t_p, t_p, d_beta);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, beta2);

  v_sub(t_gamma, alpha2, alpha);

  return Gilbert3D_xyz2d( cur_idx, q,
                          t_p,
                          t_alpha, gamma2, t_gamma );
}

//-----------------------------------------------
//    ____    __     _____     _____             
//   |_  /___/ / __ / <  / ___/ /_  |_ ____ _____
//  _/_ </ _  / / // // / / _  / __/\ \ / // /_ /
// /____/\_,_/  \___//_/  \_,_/____/_\_\\_, //__/
//                                     /___/     
//-----------------------------------------------

int Gilbert3DJ1_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma) {
  int alpha2[3],  d_alpha[3], t_alpha[3],
      beta2[3],   d_beta[3],  t_beta[3],
      gamma2[3],  d_gamma[3], t_gamma[3];

  int a,  b,  g,
      a2, b2, g2;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);
  v_div2(gamma2, gamma);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);
  g2 = abs_sum_v(gamma2);

  if ((a > 2) && ((a2 % 2) == 0)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b2 % 2) == 1)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }

  if ((g > 2) && ((g2 % 2) == 1)) {
    v_add(gamma2, gamma2, d_gamma);
    g2 = abs_sum_v(gamma2);
  }

  nxt_idx = cur_idx + (g2*a2*b2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            p,
                            gamma2, alpha2, beta2 );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (b*(g-g2)*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    v_add(t_p, p, gamma2);

    v_sub(t_beta, gamma, gamma2);
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            beta, t_beta, alpha2 );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (a*(b-b2)*g2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    v_sub(t_p, gamma2, d_gamma);
    v_add(t_p, t_p, beta);
    v_sub(t_p, t_p, d_beta);
    v_add(t_p, t_p, p);

    v_sub(t_beta, beta2, beta);

    v_neg(t_gamma, gamma2);

    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            alpha, t_beta, t_gamma );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (b*(g-g2)*(a-a2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    v_sub(t_p, alpha, d_alpha);
    v_add(t_p, t_p, beta);
    v_sub(t_p, t_p, d_beta);
    v_add(t_p, t_p, gamma2);
    v_add(t_p, t_p, p);

    v_neg(t_alpha, beta);

    v_sub(t_beta, gamma, gamma2);

    v_sub(t_gamma, alpha2, alpha);

    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            t_alpha, t_beta, t_gamma );
  }
  cur_idx = nxt_idx;

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, gamma2);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, gamma2);
  
  v_sub(t_beta, alpha2, alpha);

  return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                          t_p,
                          t_alpha, t_beta, beta2 );
}

//--------------------------------------------------
//    ____    __     _____                ___     __
//   |_  /___/ / __ / <  / __ ____ _____ |_  |___/ /
//  _/_ </ _  / / // // /  \ \ / // /_ // __// _  / 
// /____/\_,_/  \___//_/  /_\_\\_, //__/____/\_,_/  
//                            /___/                 
//--------------------------------------------------

int Gilbert3DJ1_xyz2d(int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma) {
  int alpha2[3],  d_alpha[3], t_alpha[3],
      beta2[3],   d_beta[3],  t_beta[3],
      gamma2[3],  d_gamma[3], t_gamma[3];

  int a,  b,  g,
      a2, b2, g2;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);
  v_div2(gamma2, gamma);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);
  g2 = abs_sum_v(gamma2);

  if ((a > 2) && ((a2 % 2) == 0)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b2 % 2) == 1)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }

  if ((g > 2) && ((g2 % 2) == 1)) {
    v_add(gamma2, gamma2, d_gamma);
    g2 = abs_sum_v(gamma2);
  }

  t_p[0] = p[0];
  t_p[1] = p[1];
  t_p[2] = p[2];

  if (inBounds(q, t_p, gamma2, alpha2, beta2)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            p,
                            gamma2, alpha2, beta2 );
  }
  cur_idx += (g2*a2*b2);

  v_add(t_p, p, gamma2);
  v_sub(t_beta, gamma, gamma2);
  if (inBounds(q, t_p, beta, t_beta, alpha2)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            beta, t_beta, alpha2 );
  }
  cur_idx += b*(g-g2)*a2;

  v_sub(t_p, gamma2, d_gamma);
  v_add(t_p, t_p, beta);
  v_sub(t_p, t_p, d_beta);
  v_add(t_p, t_p, p);

  v_sub(t_beta, beta2, beta);

  v_neg(t_gamma, gamma2);

  v_sub(t_beta, beta2, beta);
  v_neg(t_gamma, gamma2);
  if (inBounds(q, t_p, alpha, t_beta, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            alpha, t_beta, t_gamma );
  }
  cur_idx += a*(b-b2)*g2;

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, beta);
  v_sub(t_p, t_p, d_beta);
  v_add(t_p, t_p, gamma2);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, beta);
  v_sub(t_beta, gamma, gamma2);
  v_sub(t_gamma, alpha2, alpha);

  if (inBounds(q, t_p, t_alpha, t_beta, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            t_alpha, t_beta, t_gamma );
  }
  cur_idx += (b*(g-g2)*(a-a2));

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, gamma2);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, gamma2);
  v_sub(t_beta, alpha2, alpha);

  return Gilbert3D_xyz2d( cur_idx, q,
                          t_p,
                          t_alpha, t_beta, beta2 );
}

//-----------------------------------------------
//    ____    __     _____       _____             
//   |_  /___/ / __ / /_  |  ___/ /_  |_ ____ _____
//  _/_ </ _  / / // / __/  / _  / __/\ \ / // /_ /
// /____/\_,_/  \___/____/  \_,_/____/_\_\\_, //__/
//                                       /___/     
//-----------------------------------------------

int Gilbert3DJ2_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma) {
  int alpha2[3],  d_alpha[3], t_alpha[3],
      beta2[3],   d_beta[3],  t_beta[3],
      gamma2[3],  d_gamma[3], t_gamma[3];

  int a,  b,  g,
      a2, b2, g2;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);
  v_div2(gamma2, gamma);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);
  g2 = abs_sum_v(gamma2);

  if ((a > 2) && ((a2 % 2) == 0)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b2 % 2) == 1)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }

  if ((g > 2) && ((g2 % 2) == 1)) {
    v_add(gamma2, gamma2, d_gamma);
    g2 = abs_sum_v(gamma2);
  }

  nxt_idx = cur_idx + (b2*g*a2);
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            p,
                            beta2, gamma, alpha2 );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (g2*a*(b-b2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {
    v_add(t_p, p, beta2);
    v_sub(t_gamma, beta, beta2);
    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            gamma2, alpha, t_gamma );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (a*(b-b2)*(g-g2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    v_add(t_p, beta2, gamma2);
    v_add(t_p, t_p, p);

    v_sub(t_beta, beta, beta2);

    v_sub(t_gamma, gamma, gamma2);

    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            alpha, t_beta, t_gamma );
  }
  cur_idx = nxt_idx;

  nxt_idx = cur_idx + (b2*(g-g2)*(a-a2));
  if ((cur_idx <= dst_idx) && (dst_idx < nxt_idx)) {

    v_sub(t_p, alpha, d_alpha);
    v_add(t_p, t_p, beta2);
    v_sub(t_p, t_p, d_beta);
    v_add(t_p, t_p, gamma2);
    v_add(t_p, t_p, p);

    v_neg(t_alpha, beta2);

    v_sub(t_beta, gamma, gamma2);

    v_sub(t_gamma, alpha2, alpha);

    return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                            t_p,
                            t_alpha, t_beta, t_gamma );
  }
  cur_idx = nxt_idx;

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, gamma2);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, gamma2);

  v_sub(t_beta, alpha2, alpha);

  return Gilbert3D_d2xyz( u, dst_idx, cur_idx,
                          t_p,
                          t_alpha, t_beta, beta2);
}

//----------------------------------------------------
//    ____    __     _____                  ___     __
//   |_  /___/ / __ / /_  |  __ ____ _____ |_  |___/ /
//  _/_ </ _  / / // / __/   \ \ / // /_ // __// _  / 
// /____/\_,_/  \___/____/  /_\_\\_, //__/____/\_,_/  
//                              /___/                 
//----------------------------------------------------


int Gilbert3DJ2_xyz2d(int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma) {
  int alpha2[3],  d_alpha[3], t_alpha[3],
      beta2[3],   d_beta[3],  t_beta[3],
      gamma2[3],  d_gamma[3], t_gamma[3];

  int a,  b,  g,
      a2, b2, g2;

  int t_p[3];
  int nxt_idx;

  v_div2(alpha2, alpha);
  v_div2(beta2, beta);
  v_div2(gamma2, gamma);

  v_delta(d_alpha, alpha);
  v_delta(d_beta, beta);
  v_delta(d_gamma, gamma);

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a2 = abs_sum_v(alpha2);
  b2 = abs_sum_v(beta2);
  g2 = abs_sum_v(gamma2);

  if ((a > 2) && ((a2 % 2) == 0)) {
    v_add(alpha2, alpha2, d_alpha);
    a2 = abs_sum_v(alpha2);
  }

  if ((b > 2) && ((b2 % 2) == 1)) {
    v_add(beta2, beta2, d_beta);
    b2 = abs_sum_v(beta2);
  }

  if ((g > 2) && ((g2 % 2) == 1)) {
    v_add(gamma2, gamma2, d_gamma);
    g2 = abs_sum_v(gamma2);
  }

  t_p[0] = p[0];
  t_p[1] = p[1];
  t_p[2] = p[2];

  if (inBounds(q, p, beta2, gamma, alpha2)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            p,
                            beta2, gamma, alpha2 );
  }
  cur_idx += (b2*g*a2);

  v_add(t_p, p, beta2);
  v_sub(t_gamma, beta, beta2);
  if (inBounds(q, t_p, gamma2, alpha, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            gamma2, alpha, t_gamma );
  }
  cur_idx += (g2*a*(b-b2));

  v_add(t_p, beta2, gamma2);
  v_add(t_p, t_p, p);

  v_sub(t_beta, beta, beta2);
  v_sub(t_gamma, gamma, gamma2);
  if (inBounds(q, t_p, alpha, t_beta, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            alpha, t_beta, t_gamma );
  }
  cur_idx += (a*(b-b2)*(g-g2));

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, beta2);
  v_sub(t_p, t_p, d_beta);
  v_add(t_p, t_p, gamma2);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, beta2);
  v_sub(t_beta, gamma, gamma2);
  v_sub(t_gamma, alpha2, alpha);

  if (inBounds(q, t_p, t_alpha, t_beta, t_gamma)) {
    return Gilbert3D_xyz2d( cur_idx, q,
                            t_p,
                            t_alpha, t_beta, t_gamma );
  }
  cur_idx += (b2*(g-g2)*(a-a2));

  v_sub(t_p, alpha, d_alpha);
  v_add(t_p, t_p, gamma2);
  v_sub(t_p, t_p, d_gamma);
  v_add(t_p, t_p, p);

  v_neg(t_alpha, gamma2);

  v_sub(t_beta, alpha2, alpha);

  return Gilbert3D_xyz2d( cur_idx, q,
                          t_p,
                          t_alpha, t_beta, beta2);
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
int Gilbert3D_d2xyz(int *u, int dst_idx, int cur_idx, int *p, int *alpha, int *beta, int *gamma) {
  int a,b,g, a0,b0,g0;

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a0 = (a % 2);
  b0 = (b % 2);
  g0 = (g % 2);

  // base cases
  //
  if ((a == 2) &&
      (b == 2) &&
      (g == 2)) {
    return Hilbert2x2x2_d2xyz(u, dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if (a == 1) { return Gilbert2D_d2xyz(u, dst_idx, cur_idx, p, beta, gamma); }
  if (b == 1) { return Gilbert2D_d2xyz(u, dst_idx, cur_idx, p, alpha, gamma); }
  if (g == 1) { return Gilbert2D_d2xyz(u, dst_idx, cur_idx, p, alpha, beta); }

  // eccentric cases
  //
  if (((3*a) > (5*b)) &&
      ((3*a) > (5*g))) {
    return Gilbert3DS0_d2xyz(u, dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if (((2*b) > (3*g)) ||
      ((2*b) > (3*a))) {
    return Gilbert3DS2_d2xyz(u, dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if ((2*g) > (3*b)) {
    return Gilbert3DS1_d2xyz(u, dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  // bulk recursion
  //
  if (g0 == 0) {
    return Gilbert3DJ0_d2xyz(u, dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  if ((a0 == 0) || (b0 == 0)) {
    return Gilbert3DJ1_d2xyz(u, dst_idx, cur_idx, p, alpha, beta, gamma);
  }

  // a0 == b0 == g0 == 1
  //
  return Gilbert3DJ2_d2xyz(u, dst_idx, cur_idx, p, alpha, beta, gamma);
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
int Gilbert3D_xyz2d(int cur_idx, int *q, int *p, int *alpha, int *beta, int *gamma) {
  int a,b,g, a0,b0,g0;

  a = abs_sum_v(alpha);
  b = abs_sum_v(beta);
  g = abs_sum_v(gamma);

  a0 = (a % 2);
  b0 = (b % 2);
  g0 = (g % 2);

  // base cases
  //
  if ((a == 2) &&
      (b == 2) &&
      (g == 2)) {
    return Hilbert2x2x2_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if (a == 1) { return Gilbert2D_xyz2d(cur_idx, q, p, beta, gamma); }
  if (b == 1) { return Gilbert2D_xyz2d(cur_idx, q, p, alpha, gamma); }
  if (g == 1) { return Gilbert2D_xyz2d(cur_idx, q, p, alpha, beta); }

  // eccentric cases
  //
  if (((3*a) > (5*b)) &&
      ((3*a) > (5*g))) {
    return Gilbert3DS0_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if (((2*b) > (3*g)) ||
      ((2*b) > (3*a))) {
    return Gilbert3DS2_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if ((2*g) > (3*b)) {
    return Gilbert3DS1_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  // bulk recursion
  //
  if (g0 == 0) {
    return Gilbert3DJ0_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  if ((a0 == 0) || (b0 == 0)) {
    return Gilbert3DJ1_xyz2d(cur_idx, q, p, alpha, beta, gamma);
  }

  // a0 == b0 == g0 == 1
  //
  return Gilbert3DJ2_xyz2d(cur_idx, q, p, alpha, beta, gamma);
}

/***************************************
 *           __          __  _         
 *  ___ ____/ /__ ____  / /_(_)  _____ 
 * / _ `/ _  / _ `/ _ \/ __/ / |/ / -_)
 * \_,_/\_,_/\_,_/ .__/\__/_/|___/\__/ 
 *              /_/                    
 ***************************************/


// Gilbert3D adaptive method
//
// adapt_method is:
//
//   HARMONY(0)     : choose endpoint based on harmonious cuboid subdivision
//   HAMILTONIAN(1) : choose endpoint to make path notch-free
//   AXIS(2)        : choose endpoint on alpha axis, axis order as specified 
//
// Return:
//
//   0 on success
//
// u will hold xyz coordinates of specified index
//
int Gilbert3DAdapt_d2xyz(int *u, int idx, int width, int height, int depth, int adapt_method) {
  int w0, h0, d0, p0[3] = {0};

  int alpha[3], beta[3], gamma[3];

  w0 = (width%2);
  h0 = (height%2);
  d0 = (depth%2);

  if (adapt_method == HARMONY) {

    if ((width >= height) && (width >= depth)) {
      alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
       beta[0] = 0;      beta[1] = height;  beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }
    else if ((height >= width) && (height >= depth)) {
      alpha[0] = 0;     alpha[1] = height; alpha[2] = 0;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }
    else {
      alpha[0] = 0;     alpha[1] = 0;      alpha[2] = depth;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = height; gamma[2] = 0;
    }

  }
  
  else if (adapt_method == HAMILTONIAN) {

    if (w0 == 0) {
      alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
       beta[0] = 0;      beta[1] = height;  beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }

    else if (h0 == 0) {
      alpha[0] = 0;     alpha[1] = height; alpha[2] = 0;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }

    else if (d0 == 0) {
      alpha[0] = 0;     alpha[1] = 0;      alpha[2] = depth;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = height; gamma[2] = 0;
    }

    else {
      alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
       beta[0] = 0;      beta[1] = height;  beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }

  }

  else {

    alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
     beta[0] = 0;      beta[1] = height;  beta[2] = 0;
    gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
  }

  return Gilbert3D_d2xyz(u, idx, 0, p0, alpha, beta, gamma);
}


// Gilbert3D adaptive method
//
// q is int[3] of coordinates to convert to Gilbert index.
//
// adapt_method is:
//
//   HARMONY(0)     : choose endpoint based on harmonious cuboid subdivision
//   HAMILTONIAN(1) : choose endpoint to make path notch-free
//   AXIS(2)        : choose endpoint on alpha axis, axis order as specified 
//
// Return:
//
//   non-negative index
//
int Gilbert3DAdapt_xyz2d(int *q, int width, int height, int depth, int adapt_method) {
  int w0, h0, d0, p0[3] = {0};

  int alpha[3], beta[3], gamma[3];

  w0 = (width%2);
  h0 = (height%2);
  d0 = (depth%2);

  if (adapt_method == HARMONY) {

    if ((width >= height) && (width >= depth)) {
      alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
       beta[0] = 0;      beta[1] = height;  beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }
    else if ((height >= width) && (height >= depth)) {
      alpha[0] = 0;     alpha[1] = height; alpha[2] = 0;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }
    else {
      alpha[0] = 0;     alpha[1] = 0;      alpha[2] = depth;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = height; gamma[2] = 0;
    }

  }
  
  else if (adapt_method == HAMILTONIAN) {

    if (w0 == 0) {
      alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
       beta[0] = 0;      beta[1] = height;  beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }

    else if (h0 == 0) {
      alpha[0] = 0;     alpha[1] = height; alpha[2] = 0;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }

    else if (d0 == 0) {
      alpha[0] = 0;     alpha[1] = 0;      alpha[2] = depth;
       beta[0] = width;  beta[1] = 0;       beta[2] = 0;
      gamma[0] = 0;     gamma[1] = height; gamma[2] = 0;
    }

    else {
      alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
       beta[0] = 0;      beta[1] = height;  beta[2] = 0;
      gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
    }

  }

  else {
    alpha[0] = width; alpha[1] = 0;      alpha[2] = 0;
     beta[0] = 0;      beta[1] = height;  beta[2] = 0;
    gamma[0] = 0;     gamma[1] = 0;      gamma[2] = depth;
  }

  return Gilbert3D_xyz2d(0, q, p0, alpha, beta, gamma);
}

/***************************************
 *           __          __  _         
 *  ___ ____/ /__ ____  / /_(_)  _____ 
 * / _ `/ _  / _ `/ _ \/ __/ / |/ / -_)
 * \_,_/\_,_/\_,_/ .__/\__/_/|___/\__/ 
 *              /_/                    
 ***************************************/



/***************************************************
 *   ______ ____           __  ____ ___    __    __ 
 *  / ___(_) / /  ___ ____/ /_|_  // _ \__/ /___/ /_
 * / (_ / / / _ \/ -_) __/ __//_ </ // /_  __/_  __/
 * \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
 *                                                  
 ***************************************************/

#define GILBERTPP_MAIN
#ifdef GILBERTPP_MAIN

#include <string.h>

/*********************
 *               _
 *   __ _  ___ _(_)__
 *  /  ' \/ _ `/ / _ \
 * /_/_/_/\_,_/_/_//_/
 *
 **********************/


int main(int argc, char **argv) {
  int w, h, d;
  int x, y, z;
  int idx;

  int xyz[3] = {0}, p[3] = {0},
      alpha[3], beta[3], gamma[3];

  char buf[1024];

  int adapt_method = HARMONY;
  int r;

  w = 1;
  h = 1;
  d = 1;

  if (argc < 4) {
    printf("provide args\n");
    printf("\n");
    printf("usage:\n");
    printf("\n");
    printf("  gilbert3dpp <op> <width> <height> [depth]\n");
    printf("\n");
    printf("    op      - one of \"xy2d\",\"2dxy\",\"xyz2d\",\"d2xyz\"\n");
    printf("    depth   - default to 1 for 3D Gilbert++ with no depth specified\n");
    printf("\n");
    exit(-1);
  }

  strncpy(buf, argv[1], 1023);
  buf[1023]='\0';

  w = atoi(argv[2]);
  h = atoi(argv[3]);
  if (argc > 4) {
    d = atoi(argv[4]);
  }

  if ((w <= 0) || (h <= 0) || (d <= 0)) {
    fprintf(stderr, "width, height, depth must all be greater than 0\n");
    exit(-1);
  }

  if      (strncmp("d2xyz.0", buf, 1023) == 0) { adapt_method = HARMONY; }
  else if (strncmp("d2xyz.1", buf, 1023) == 0) { adapt_method = HAMILTONIAN; }
  else if (strncmp("d2xyz.2", buf, 1023) == 0) { adapt_method = AXIS; }

  else if (strncmp("xyz2d.0", buf, 1023) == 0) { adapt_method = HARMONY; }
  else if (strncmp("xyz2d.1", buf, 1023) == 0) { adapt_method = HAMILTONIAN; }
  else if (strncmp("xyz2d.2", buf, 1023) == 0) { adapt_method = AXIS; }

  else if (strncmp("d2xy.0", buf, 1023) == 0) { adapt_method = HARMONY; }
  else if (strncmp("d2xy.1", buf, 1023) == 0) { adapt_method = HAMILTONIAN; }
  else if (strncmp("d2xy.2", buf, 1023) == 0) { adapt_method = AXIS; }

  else if (strncmp("xy2d.0", buf, 1023) == 0) { adapt_method = HARMONY; }
  else if (strncmp("xy2d.1", buf, 1023) == 0) { adapt_method = HAMILTONIAN; }
  else if (strncmp("xy2d.2", buf, 1023) == 0) { adapt_method = AXIS; }

  p[0] = 0;
  p[1] = 0;
  p[2] = 0;

  if (strncmp("xy2d", buf, 1023) == 0) {

    alpha[0] = w; alpha[1] = 0; alpha[2] = 0;
     beta[0] = 0;  beta[1] = h;  beta[2] = 0;

    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = 0;
        idx = Gilbert2D_xyz2d( 0, xyz, p, alpha, beta );
        printf("%i %i %i\n", idx, x, y);
      }
    }

  }
  else if (strncmp("d2xy", buf, 1023) == 0) {

    alpha[0] = w; alpha[1] = 0; alpha[2] = 0;
     beta[0] = 0;  beta[1] = h;  beta[2] = 0;

    for (idx = 0; idx < (w*h); idx++) {
      Gilbert2D_d2xyz( xyz, idx, 0, p, alpha, beta );
      printf("%i %i\n", xyz[0], xyz[1]);
    }

  }
  else if (strncmp("xyz2d", buf, 1023) == 0) {

    alpha[0] = w; alpha[1] = 0; alpha[2] = 0;
     beta[0] = 0;  beta[1] = h;  beta[2] = 0;
    gamma[0] = 0; gamma[1] = 0; gamma[2] = d;

    for (z = 0; z < d; z++) {
      for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {

          xyz[0] = x;
          xyz[1] = y;
          xyz[2] = z;

          idx = Gilbert3D_xyz2d( 0, xyz, p, alpha, beta, gamma );
          printf("%i %i %i %i\n", idx, x, y, z);
        }
      }
    }

  }

  else if ((strncmp("xyz2d.0", buf, 1023) == 0) ||
           (strncmp("xyz2d.1", buf, 1023) == 0) ||
           (strncmp("xyz2d.2", buf, 1023) == 0)) {

    for (z = 0; z < d; z++) {
      for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
          xyz[0] = x;
          xyz[1] = y;
          xyz[2] = z;
          idx = Gilbert3DAdapt_xyz2d( xyz, w, h, d, adapt_method );
          printf("%i %i %i %i\n", idx, x, y, z);
        }
      }
    }

  }


  else if (strncmp("d2xyz", buf, 1023) == 0) {

    alpha[0] = w; alpha[1] = 0; alpha[2] = 0;
     beta[0] = 0;  beta[1] = h;  beta[2] = 0;
    gamma[0] = 0; gamma[1] = 0; gamma[2] = d;

    for (idx = 0; idx < (w*h*d); idx++) {
      r = Gilbert3D_d2xyz( xyz, idx, 0, p, alpha, beta, gamma );
      printf("%i %i %i\n", xyz[0], xyz[1], xyz[2]);
    }

  }

  else if ((strncmp("d2xyz.0", buf, 1023) == 0) ||
           (strncmp("d2xyz.1", buf, 1023) == 0) ||
           (strncmp("d2xyz.2", buf, 1023) == 0)) {

    for (idx = 0; idx < (w*h*d); idx++) {
      r = Gilbert3DAdapt_d2xyz( xyz, idx, w, h, d, adapt_method );
      printf("%i %i %i\n", xyz[0], xyz[1], xyz[2]);
    }

  }


  exit(0);
}

#endif



