Gilbert Extensions
===

This directory holds reference implementations for the extensions to the original Gilbert curve,
called `Gilbert++`.
These extensions are an accompaniment to the Gilbert curve paper (see [jakubcerveny/gilbet-paper](https://github.com/jakubcerveny/gilbert-paper/)
for details).
Please refer to the paper for details on the underlying algorithm.

There are reference implementations for three languages, `C`, `Python` and `JavaScript`.

| Program | Function | Type | Output |
|---------|----------|------|--------|
| `gilbert3dpp.py` | `Gilbert2D(width,height)` | Asynchronous | `(x,y,z)` |
| `gilbert3dpp.py` | `Gilbert2DAdapt(width,height)` | Asynchronous | `(x,y,z)` |
| `gilbert3dpp.py` | `Gilbert3D(width,height,depth)` | Asynchronous | `(x,y,z)` |
| `gilbert3dpp.py` | `Gilbert3DAdapt(width,height,depth)` | Asynchronous | `(x,y,z)` |
| `gilbert3dpp.py` | `Gilbert2D_xyz2d(idx,xyz,p0,alpha,beta)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.py` | `Gilbert2D_d2xyz(idx,startIndex,p0,alpha,beta)` | Synchronous, Random Access | `(x,y,z)` |
| `gilbert3dpp.py` | `Gilbert3D_xyz2d(idx,xyz,p0,alpha,beta,gamma)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.py` | `Gilbert3D_d2xyz(idx,startIndex,p0,alpha,beta,gamma)` | Synchronous, Random Access | `(x,y,z)` |
| `gilbert3dpp.js` | `Gilbert2DAsync(width,height)` | Asynchronous | `(x,y,z)` (generator) |
| `gilbert3dpp.js` | `Gilbert3DAsync(width,height,depth)` | Asynchronous | `(x,y,z)` (generator) |
| `gilbert3dpp.js` | `Gilbert2D_xyz2d(startIndex,xyz,alpha,beta,gamma)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.js` | `Gilbert2DAdapt_xy2d(xy,width,height)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.js` | `Gilbert2D_d2xyz(idx,startIndex,p0,alpha,beta,gamma)` | Synchronous, Random Access | `(x,y,z)` |
| `gilbert3dpp.js` | `Gilbert3D_xyz2d(startIndex,xyz,alpha,beta,gamma)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.js` | `Gilbert3DAdapt_xy2d(xy,width,height)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.js` | `Gilbert3D_d2xyz(idx,startIndex,p0,alpha,beta,gamma)` | Synchronous, Random Access | `(x,y,z)` |
| `gilbert3dpp.c` | `Gilbert2D_xyz2d(idx,xyz,p0,alpha,beta)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.c` | `Gilbert2D_d2xyz(xyz,idx,startIndex,p0,alpha,beta)` | Synchronous, Random Access | `(x,y,z)` (in `xyz`) |
| `gilbert3dpp.c` | `Gilbert3D_xyz2d(idx,xyz,p0,alpha,beta,gamma)` | Synchronous, Random Access | `index` |
| `gilbert3dpp.c` | `Gilbert3D_d2xyz(xyz,idx,startIndex,p0,alpha,beta,gamma)` | Synchronous, Random Access | `(x,y,z)` (in `xyz`) |


Each program can be run from the command line.
Source code prioritizes ease of reading over efficiency.

Adaptive Methods
---

The adaptive methods have an option to specify what strategy should be employed when choosing an endpoint.

There are three methods:

| Name | Integer Code | Description |
|------|--------------|----|
| Harmony | 0 | Prioritizes a "harmonious" subdivision and chooses the endpoint on the axis that is initially largest |
| Hamiltonian | 1 | Prioritizes a notch-free Hamiltonian path, choosing from the three axis-aligned neighboring endpoints, from the starting point, on the cuboid edge |
| Axis | 2 | Uses axis order specified (default for other, non-adaptive, functions above) |


Tests
---

Tests can be run in the `tests/` directory.
From here:

```
$ cd ../tests/
$ ./runtests-extensions.sh
```

