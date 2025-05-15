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

/* ****
   __       __               ___              __  _
  / /  ___ / /__  ___ ____  / _/_ _____  ____/ /_(_)__  ___  ___
 / _ \/ -_) / _ \/ -_) __/ / _/ // / _ \/ __/ __/ / _ \/ _ \(_-<
/_//_/\__/_/ .__/\__/_/   /_/ \_,_/_//_/\__/\__/_/\___/_//_/___/
          /_/
**** */


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
    d = alpha[0] + beta[0] + gamma[0];

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

/* ****
   __       __               ___              __  _
  / /  ___ / /__  ___ ____  / _/_ _____  ____/ /_(_)__  ___  ___
 / _ \/ -_) / _ \/ -_) __/ / _/ // / _ \/ __/ __/ / _ \/ _ \(_-<
/_//_/\__/_/ .__/\__/_/   /_/ \_,_/_//_/\__/\__/_/\___/_//_/___/
          /_/
**** */

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


#define GILBERTPP_MAIN
#ifdef GILBERTPP_MAIN

#include <string.h>

/* ***
              _
  __ _  ___ _(_)__
 /  ' \/ _ `/ / _ \
/_/_/_/\_,_/_/_//_/

*** */


int main(int argc, char **argv) {
  int w, h, d;
  int x, y, z;
  int idx;

  int xyz[3] = {0}, p[3] = {0},
      alpha[3], beta[3], gamma[3];

  char buf[1024];

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
    exit(-1);
  }

  if (strncmp("xy2d", buf, 1023) == 0) {

    p[0] = 0;
    p[2] = 0;
    p[2] = 0;

    alpha[0] = w;
    alpha[1] = 0;
    alpha[2] = 0;

    beta[0] = 0;
    beta[1] = h;
    beta[2] = 0;

    for (x = 0; x < w; x++) {
      for (y = 0; y < h; y++) {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = 0;
        idx = Gilbert_xy2d( x, y, w, h );
        printf("%i %i %i\n", idx, x, y);
      }
    }

  }
  else if (strncmp("d2xy", buf, 1023) == 0) {

    p[0] = 0;
    p[2] = 0;
    p[2] = 0;

    alpha[0] = w;
    alpha[1] = 0;
    alpha[2] = 0;

    beta[0] = 0;
    beta[1] = h;
    beta[2] = 0;

    for (idx = 0; idx < (w*h); idx++) {
      Gilbert2D_d2xyz( xyz, idx, 0, p, alpha, beta );
      printf("%i %i\n", xyz[0], xyz[1]);
    }

  }
  else if (strncmp("xyz2d", buf, 1023) == 0) {

    for (x = 0; x < w; x++) {
      for (y = 0; y < h; y++) {
        for (z = 0; z < d; z++) {
          //idx = Gilbert_xyz2d( x,y,z, w,h,d );
          //printf("%i %i %i %i\n", idx, x, y, z);
        }
      }
    }

  }

  else if (strncmp("d2xyz", buf, 1023) == 0) {



    for (idx = 0; idx < (w*h*d); idx++) {
      //Gilbert_d2xyz( xyz, idx, 0, p, alpha, beta )
      //printf("%i %i %i\n", x, y, z);
    }

  }


  exit(0);
}

#endif
