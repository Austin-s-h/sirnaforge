"""Fixed-length benchmark artifacts (#109).

A benchmark artifact is a versioned, on-disk directory pairing a measured siRNA/miRNA efficacy
panel with the current fixed-length design/filter path, without adding a new execution path or
touching any global default. ``sirnaforge.benchmark.artifact`` owns the artifact/manifest schema
and the lossless CSV/JSON round trip; later slices (panel registry, ``prepare``/``design``
commands) add their own modules under this package.

Deliberately no re-exports here. A symbol added to this docstring's module and not to an import
statement below it is a promise this package does not keep -- so there is no import statement
below it: every consumer names the full module path (``sirnaforge.benchmark.artifact.X``), and a
slice that only owns one file in this package cannot be broken by another slice's edit to this one.
"""
