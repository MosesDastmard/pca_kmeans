
function  B=NNW(X,XM)

[n,J]=size(X);
vecXM=XM(:);
vecX=X(:);
b= pinv(vecXM)*vecX;
B=diag(b);
%%%%%%%%%%%%%%%%%%%%%%%%%%


  