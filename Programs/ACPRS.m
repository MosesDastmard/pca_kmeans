%%%%%%%%%%%%%%%%%%
% ACP Restricted %
%%%%%%%%%%%%%%%%%%

% Maurizio Vichi 
% October 2012
%
function [V, Ars, Yrs]=ACPRS(X);
% initialize variables
[n,J]=size(X);

uno=ones(J,1);
C=eye(J);
% centring matrix
%Jm=eye(n)-(1/n)*ones(n);
Xs=zscore(X,1);
% compute var-covar matrix 

%S=(1/n)*X'*Jm*X;
% Standardize data
R=cov(Xs,1);

%Xs=Jm*X*diag(diag(S))^-0.5;
V=[];
% compute var-covar matrix of standardized data
%R=(1/n)*Xs'*Xs;
% compute ACP
[A,C]=eig(R);
[c,ic]=sort(diag(C), 'descend');
C=diag(c);
A=A(:,ic);
inc=find(diag(C)>=1);
AA=C.^0.5*A(:,1:size(inc,1));
[vm, m]=max(abs(AA'));
Ars=zeros(J,size(inc,1));
V=zeros(J,size(inc,1));
for j=1:J
    Ars(j,m(j))=AA(j,m(j));
    V(j,m(j))=1;
end
%Ars=Ars*(diag(sum(Ars.^2))^-0.5);
Ars=Ars*inv(Ars'*Ars)^0.5;
Yrs=Xs*Ars;
