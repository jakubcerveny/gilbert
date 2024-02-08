

# gilbert

Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of
arbitrary (non-power of two) sizes.

The discrete [Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve) is a
widely used space-filling curve to map between N-dimensional and 1-D spaces
while preserving locality. However, classical algorithms only work for domains
whose sides are powers of two.

We present a simple recursive algorithm that generalizes the Hilbert curve
to rectangles of arbitrary sizes in 2D, and cuboids of even sizes in 3D.

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/55x31.png)


### Previous Work

A 2D algorithm that combines Peano (3x3) and Hilbert (2x2) blocks was published
by Lutz Tautenhahn in 2003:

> [1] Lutz Tautenhahn: Draw a Space-Filling Curve of Arbitrary Size, http://lutanho.net/pic2html/draw_sfc.html, 2003.

However, this algorithm is complex and would be very difficult to generalize to
3D. Another approach that focuses on evaluation speed is described in

> [2] Zhang J., Kamata S., Ueshige Y.: A Pseudo-Hilbert Scan Algorithm for Arbitrarily-Sized Rectangle Region, IWICPAS 2006.

This method generalizes to 3D, but it is again complex and the resulting curves
are a bit uneven.


### Our Algorithm

The idea is to recursively apply the following template to obtain a Hilbert-like
curve. A general rectangle with a known orientation is split into three regions
("up", "right", "down"), for which the function calls itself recursively, until
a trivial path can be produced.

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/algorithm.svg?sanitize=true)

The trick is to avoid odd sizes when halving the dimensions (a/2, b/2). The
reason is that if the width of a rectangle is odd and its height is even, it is
impossible to generate a continuous path that goes from the lower left corner to
the lower right corner. In such a case, at least one "diagonal" step needs to be
inserted, and the algorithm tries to avoid that whenever it can. Since there is
usually some freedom in positioning the horizontal split, the algorithm
increments the length "b/2" by one if it happens to be odd. That way, the "up"
and "down" rectangles are easier to fill correctly.

There is also the case of a long "right" step, which is resolved by two
recursions only in the same direction. In this case, the length "a/2" is chosen
to be even. If there ever is a situation where the diagonal step cannot be
avoided, it will happen in the top right corner of the original domain.


### Examples

Running `gilbert2d` with two arguments (width, height) produces a space-filling
curve with orthogonal steps only, as long as the width is even (100x63 shown):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/100x63.png)

If the sizes are powers of two, a standard Hilbert curve is generated.
The algorithm extends naturally to 3D (8x6x4):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/8x6x4.png)

40x30x20:

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/40x30x20.png)

Very flat is OK too (20x12x2):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/20x12x2.png)


### Odd Sizes

In 2D, if the larger dimension is odd and the smaller is even, a single diagonal
step cannot be avoided (15x12):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/15x12.png)

In 3D this is much worse, so odd dimensions should not be used (7x6x4):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/7x6x4.png)


### Visualizing the Results

A simple Octave script (`plotpath.m`) is included to help visualize the
coordinates produced by `gilbert2d` and `gilbert3d`. The above figures were
obtained with

```
./gilbert2d.py 100 63 | octave --eval 'waitfor(plotpath(dlmread(stdin())));'
```

### Random Access Functions

The following programs provide a reference implementation to efficiently convert
from an index point along the 1D generalized Hilbert curve to spatial coordinates, `d2xy`/`d2xyz`,
and from spatial coordinates to the position along the generalized Hilbert curve, `xy2d`/`xyz2d`:

| Program | Function | Output |
|---|---|---|
| `gilbert_d2xy.py` | `gilbert_d2xy(index,width,height)` | `(x,y)` |
| `gilbert_xy2d.py` | `gilbert_xy2d(x,y,width,height)` | `index` |
| `gilbert_d2xyz.py` | `gilbert_d2xyz(index,width,height,depth)` | `(x,y,z)` |
| `gilbert_xyz2d.py` | `gilbert_xyz2d(x,y,z,width,height,depth)` | `index` |

Each of the programs run on the command line enumerate through the points or indices for their
respective `d2xy[z]` or `xy[z]2d` functions.

The `gilbert_d2xy.py` and `gilbert_d2xyz.py` programs will print out the `index`, `x`, `y` and `z` point,
where appropriate, in `x`, `y`, `z` order.

For example:

```
./gilbert_xyz2d.py 20 12 2
0 0 0 0
1 0 0 1
7 0 1 0
6 0 1 1
8 0 2 0
...
```

To recover the index order traversal, one can sort on the first element (for example `./gilbert_xyz2d.py 20 12 2 | sort -n | cut -f2- -d' '`).

Each of these functions follows the reference `gilbert2d.py` and `gilbert3d.py` implementations but
"short circuit" the recursion when the index or spatial point wouldn't land in the sub problem.
Each leg of the recursion is restricted to another rectangle or cuboid region, allowing for an easy determination
of a spatial bounds check and to count the length of the path in the sub region.
This provides an $O( \lg(N = W \cdot H \cdot D) )$ algorithm.

### Ports

There are reference implementations for JavaScript and C in the `ports` directory.

| Language | Program | Function | Output |
|----------|---------|----------|--------|
| `C`      | `gilbert.c` | `int gilbert_d2xy(int *x, int *y, int idx, int width, int height)` | Values stored in `x`,`y` |
| `C`      | `gilbert.c` | `int gilbert_xy2d(int x, int y, int width, int height)` | Returns index |
| `C`      | `gilbert.c` | `int gilbert_d2xyz(int *x, int *y int *z, int idx, int width, int height, int depth)` | Values stored in `x`,`y`,`z` |
| `C`      | `gilbert.c` | `int gilbert_xyz2d(int x, int y, int z, int width, int height, int depth)` | Returns index |
| `JS`     | `gilbert.js` | `gilbert.d2xy(idx,width,height)` | Returns `{"x":x,"y":y}` object |
| `JS`     | `gilbert.js` | `gilbert.xy2d(x,y,width,height)` | Returns index |
| `JS`     | `gilbert.js` | `gilbert.d2xyz(idx,width,height,depth)` | Returns `{"x":x,"y":y, "z":z}` object |
| `JS`     | `gilbert.js` | `gilbert.xyz2d(x,y,z,width,height,depth)` | Returns index |

### [Web Demo](https://jakubcerveny.github.io/gilbert/demo)

There is a web demo program, in the `demo/` sub directory, that can be used to to interactively explore
2D generalized Hilbert curves.

A live version can be found [here](https://jakubcerveny.github.io/gilbert/demo).

---

Author: Jakub Červený. This [code](https://github.com/jakubcerveny/gilbert) is released under the 2-clause BSD license.

