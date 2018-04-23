# PQSQ-DataApproximators

This is a set of Matlab procedures for performing PQSQ-based data approximation.

PQSQ stands for "piece-wise quadratic sub-quadratic" error function which can approximate a large family of error functions in any standard machine learning algorithm and substitute the standard quadratic error function. This is a way to construct very fast and relatively accurate approximators or regressions with non-quadratic error function or with non-quadratic regularizers.

[The theory behind PQSQ methods](http://www.math.le.ac.uk/people/ag153/homepage/GorbanMirkesZinovyevNN2016.pdf).

In particular,

PQSQmean - computes the mean value with PQSQ approximation error 

pcaPQSQ - computes PCA with PQSQ approximation error 

kmeansPQSQ - computes k-means clustering using PQSQ potential

Simplest examples of use are provided in the comments to the corresponding functions

'test_data' folder contains real and synthetic data and the code used to benchmark PQSQ algorithms 
against existing L1-based PCA implementations
