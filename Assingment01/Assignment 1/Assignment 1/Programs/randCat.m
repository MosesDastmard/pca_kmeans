function [x]=mvmnCat(n,Nj,p)

% generates a random partition of n objects in Nj classes
%
% n = number of objects
% Nj = number of categories
% p = vector of Nj probabilities
%
x_=[1:Nj]';
U = multrnd(1,p,n);
x = U*x_;
