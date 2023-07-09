# TOP2

This is a script quantify the Hi-C interaction change locally between samples.

First, distance-normalized matrices (observed/expected matrices) were calculated by dividing each matrix by the average of all entries at the same distance. Then, each log2 fold change w×w submatrix was calculated at each bin along the diagonal of the whole chromosome matrix. Last, the differential score was defined by using mean of log2 submatrix.
