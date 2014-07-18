mars
====

mars is a graph layout tool for large graph visualization. By constructing a low-rank approximation to the all-pairs shortest path matrx, mars can compute stress majorization based layouts on graphs with as many as 400,000 nodes. This is well beyond the limits of standard stress majorziation layout algorithms, including those implemented in neato. It is based on the work in [Drawing Large Graphs by Low-Rank Stress Majorization](http://www.cs.berkeley.edu/~khoury/mars.pdf).

![finance256](./finance256.gif)

Build
=====

mars is dependent upon cgraph, which is part of the Graphviz library, and LAPACK. On a Mac these libraries can be installed by running
