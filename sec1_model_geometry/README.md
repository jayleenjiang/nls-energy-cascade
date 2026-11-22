# Fig. 1: resonant-family embedding (Section 3, introduction)

- `geometry/main.cpp` builds the first five generations of the CKSTT resonant set with exact rational arithmetic (Boost.Rational) and writes, for each generation step, the parents and children of every nuclear family: `geometry/rectangles_N5_R100_f<k>_to_f<k+1>.csv`. Build with `geometry/CMakeLists.txt`.
- `geometry/verify_rectangles.py` checks that every family is a rectangle (the resonance conditions); `geometry/plot_rectangles.py` is a quick plot.
- `figures/make_figure.py [OUT]` draws `cascade_embedding.pdf/.png` from the four CSV files. `figures/cascade_embedding_paper.pdf` is the file used in the paper (the paper includes a PNG export of it); the redrawn figure has the same content and a slightly different page size.
