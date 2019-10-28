

# gilbert

Generalized Hilbert ("gilbert") discrete space-filling curve for rectangular
domains of arbitrary (non-power of two) sizes.

The discrete [Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve) is a
widely used space-filling curve to map between N-dimensional and 1-D spaces
while preserving locality. However, the classical algorithms only work for
domains whose sides are powers of two.

We present a simple recursive algorithm that generalizes the Hilbert curve
to rectangles of arbitrary sizes in 2D, and cuboids of even sizes in 3D.

![55x31](https://raw.githubusercontent.com/jakubcerveny/gilbert/img/55x31.png)


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



```
./gilbert2d 15 12 | octave --eval 'waitfor(plotpath(dlmread(stdin())));'
```


TODO

