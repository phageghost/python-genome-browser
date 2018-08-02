1. Ensure compatibility with all versions of Ensembl GFF3
2. Add GTF support
3. Query overlaps using intervaloverlaps and remove requirement for intervaltrees
4. Improve speed and memory requirements by either:
    1. Using pandas to read GFF3 file instead of pure python parsing
    2. Index GFF3 file on disk and retrieve on the fly
        * Optionally use BGZIP indexing
    3. Binary search file on disk like we do for tag directories
5. Add support for stranded data
6. Add flag for fill-beneath and make default.
7. Shrink logo.
8. Look into eliminating requirement for seaborn (at this point I think we only use the styles).
    