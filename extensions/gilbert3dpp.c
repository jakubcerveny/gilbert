/*

 To the extent possible under law, the person who associated CC0 with
 this project has waived all copyright and related or neighboring rights
 to this project.

 You should have received a copy of the CC0 legalcode along with this
 work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

*/

#include <stdio.h>
#include <stdlib.h>

enum GILBERT3DPP_ADAPT_METHOD = {
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

int inBounds(int *q, int *p, int *alpha, int *beta, int *_gamma=NULL) {
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



/* ****
   __       __               ___              __  _
  / /  ___ / /__  ___ ____  / _/_ _____  ____/ /_(_)__  ___  ___
 / _ \/ -_) / _ \/ -_) __/ / _/ // / _ \/ __/ __/ / _ \/ _ \(_-<
/_//_/\__/_/ .__/\__/_/   /_/ \_,_/_//_/\__/\__/_/\___/_//_/___/
          /_/
**** */


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

    for (x = 0; x < w; x++) {
      for (y = 0; y < h; y++) {
        idx = Gilbert_xy2d( x, y, w, h );
        printf("%i %i %i\n", idx, x, y);
      }
    }

  }
  else if (strncmp("d2xy", buf, 1023) == 0) {

    for (idx = 0; idx < (w*h); idx++) {
      Gilbert_d2xy( &x, &y, idx, w, h );
      printf("%i %i\n", x, y);
    }

  }
  else if (strncmp("xyz2d", buf, 1023) == 0) {

    for (x = 0; x < w; x++) {
      for (y = 0; y < h; y++) {
        for (z = 0; z < d; z++) {
          idx = Gilbert_xyz2d( x,y,z, w,h,d );
          printf("%i %i %i %i\n", idx, x, y, z);
        }
      }
    }

  }

  else if (strncmp("d2xyz", buf, 1023) == 0) {

    for (idx = 0; idx < (w*h*d); idx++) {
      Gilbert_d2xyz( &x,&y,&z, idx, w,h,d );
      printf("%i %i %i\n", x, y, z);
    }

  }


  exit(0);
}

#endif
