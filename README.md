

# gilbert

Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of
arbitrary (non-power of two) sizes.

The discrete [Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve) is a
widely used space-filling curve to map between N-dimensional and 1-D spaces
while preserving locality. However, the classical algorithms only work for
domains whose sides are powers of two.

We present a simple recursive algorithm that generalizes the Hilbert curve
to rectangles of arbitrary sizes in 2D, and cuboids of even sizes in 3D.

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/55x31.png)


### Previous Work

A 2D algorithm that combines Peano (3x3) and Hilbert (2x2) blocks was published
by Lutz Tautenhahn in 2003:

> [1] Lutz Tautenhahn: Draw a Space-Filling Curve of Arbitrary Size, http://lutanho.net/pic2html/draw_sfc.html, 2003.

However, this algorithm is complex and would be very difficult to generalize to
3D. Another aproach that focuses on evaluation speed is described in

> [2] Zhang J., Kamata S., Ueshige Y.: A Pseudo-Hilbert Scan Algorithm for Arbitrarily-Sized Rectangle Region, IWICPAS 2006.

This method generalizes to 3D, but it is again complex and the resulting curves
are a bit uneven.


### Our Algorithm

The idea is to recursively apply the following template to obtain a Hilbert-like
curve. The trick is to avoid odd sizes when halving the dimensions (a/2, b/2).
The horizontal dimension of the template on the left should be even, if possible.

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/algorithm.svg?sanitize=true)

Running `gilbert2d` with two arguments (width, height) produces a space-filling
curve with orthogonal steps only, as long as the width is even (100x63 shown):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/100x63.png)

The algorithm extends naturally to 3D (8x6x4):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/8x6x4.png)

40x30x20:

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/40x30x20.png)

Very flat is OK too (20x12x2):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/20x12x2.png)


### Odd Sizes

In 2D, if the larger dimension is odd and the smaller is even, a single diagonal
step needs to be inserted:

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/15x12.png)

In 3D this is much worse, so odd dimensions should be avoided altogether (7x6x4):

![](https://raw.githubusercontent.com/jakubcerveny/gilbert/master/img/7x6x4.png)


### Visualizing the Results

A simple Octave script is included to help visualize the coordinates produced by
`gilbert2d` and `gilbert3d`. The above figures were obtained with this command
line:

```
./gilbert2d 100 63 | octave --eval 'waitfor(plotpath(dlmread(stdin())));'
```



---

TODO

