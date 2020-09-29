% calcola una cluster gerarchica
%       'single'    --- nearest distance
%       'complete'  --- furthest distance
%       'average'   --- unweighted average distance (UPGMA) (also known as
%                       group average)
%       'weighted'  --- weighted average distance (WPGMA)
%       'centroid'  --- unweighted center of mass distance (UPGMC) (*)
%       'median'    --- weighted center of mass distance (WPGMC) (*)
%       'ward'      --- inner squared distance (min variance algorithm) (*)
 
%    (*) The output Z for the centroid, median, and Ward's methods is
%    meaningful only if Y contains Euclidean distances.
 
%   Y = PDIST(X, DISTANCE) computes Y using DISTANCE.  Choices are:
 
%        'euclidean'   - Euclidean distance
%        'seuclidean'  - Standardized Euclidean distance, each coordinate
%                        in the sum of squares is inverse weighted by the
%                        sample variance of that coordinate
%        'cityblock'   - City Block distance
%        'mahalanobis' - Mahalanobis distance
%        'minkowski'   - Minkowski distance with exponent 2
%        'cosine'      - One minus the cosine of the included angle
%                        between observations (treated as vectors)
%        'correlation' - One minus the sample linear correlation between
%                        observations (treated as sequences of values).
%        'spearman'    - One minus the sample Spearman's rank correlation
%                        between observations (treated as sequences of values).
%        'hamming'     - Hamming distance, percentage of coordinates
%                        that differ
%        'jaccard'     - One minus the Jaccard coefficient, the
%                        percentage of nonzero coordinates that differ
%        'chebychev'   - Chebychev distance (maximum coordinate difference)
%        function      - A distance function specified using @, for
%                        example @DISTFUN

function [H,T]=gerarchica(X)
D= PDIST(X, 'euclidean');
Z= linkage(D,'average'); 
[H, T] = dendrogram(Z, 'colorthreshold',2);

