%%%%%%%%%%%%%%%%%
% ACP Step Wise %
%%%%%%%%%%%%%%%%%

% Maurizio Vichi 
% October 2012
%
function [V, Asw, Ysw]=ACPSW(X);
% initialize variables
[n,J]=size(X);
B=eye(J);
uno=ones(J,1);
C=eye(J);
% centring matrix
Jm=eye(n)-(1/n)*ones(n);

% compute var-covar matrix 
S=(1/n)*X'*Jm*X;
% Standardize data

Xs=Jm*X*diag(diag(S))^-0.5;
V=[];
% compute var-covar matrix of standardized data
R=(1/n)*Xs'*Xs;
for k=1:J;
    for j=sum(diag(B)):-1:1
        % compute ACP
        R1=B*R*B;
        [A,C]=eig(R1);
        [c,ic]=sort(diag(C), 'descend');
        C=diag(c);
        A=A(:,ic);
        if C(2,2)<1
            
            break
        end
        inc=find(diag(C)>=1);
        [vm, m]=max(abs(A(:,size(inc,1))));
        B(m,m)=0;    
    end
V=[V diag(B)];
R1=B*R*B;
[A,C]=eig(R1);
[c,ic]=sort(diag(C), 'descend');
C=diag(c);
A=A(:,ic);
Asw(:,k)=A(:,1);
Ysw(:,k)=Xs*Asw(:,k);
B=diag((uno-diag(B)));
uno=uno-V(:,k);
if sum(uno) <= 1
    if sum(uno) == 1
        V(:,k+1)=uno;
        Asw(:,k+1)=uno;
        Ysw(:,k+1)=Xs*Asw(:,k+1);
    end
    break
end
end