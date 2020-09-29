function [U,X,Xmean]=mcgen(nv,dc,nvn);
%
% random generation of cluster structures (Milligan & Cooper, JoC 1985)
%
% 21/01/2001
%
% n  = number of units
% nk = number of clusters
%
n=sum(dc);
nk=max(size(dc));
%
Xmean=zeros(nk,nv);
sd=zeros(nk,nv);
for m=1:nv,   
   % generate boundaries
   len=10+(30/m)*rand(nk,1);
   sd(:,m)=len./3;              
   Xmean(1,m)=len(1)./2;
   for k=2:nk
      Xmean(k,m)=Xmean(k-1,m)+(len(k-1)+(.5+rand(1))*(sd(k-1)+sd(k))+len(k))./2;
   end
   p=randperm(nk);        % permute centroids and sd
   sd(:,m)=sd(p,m);
   Xmean(:,m)=Xmean(p,m); 
end

X=zeros(n,nv);
U=zeros(n,nk);
h=0;
for k=1:nk,
   for i=1:dc(k)
      h=h+1;
      for j=1:nv,
         rn=2;
         while abs(rn) > 1.5
            rn=randn(1);
         end
         X(h,j)=sd(k,j)*rn;
      end
      U(h,k)=1;
   end
end
X=preproa(U*Xmean+X,1);
switch nvn,
    case 1,
    case 2,
        ind=[7 8];
        ind=[ind, ind+8, ind+16, (1:8)+24];
        X(:,ind)=X(:,ind)-U*pinv(U'*U)*U'*X(:,ind);
        X=preproa(X,1);
    case 3,
        ind=[5 6 7 8];
        ind=[ind, ind+8, ind+16, (1:8)+24];
        X(:,ind)=X(:,ind)-U*pinv(U'*U)*U'*X(:,ind);
        X=preproa(X,1);
    case 4,
        ind=[3 4 5 6 7 8];
        ind=[ind, ind+8, ind+16, (1:8)+24];
        X(:,ind)=X(:,ind)-U*pinv(U'*U)*U'*X(:,ind);
        X=preproa(X,1);    
    case 5,
        E=preproa(randn(n,32),1);
        rob=[4.36 4.36 4.36 4.36 4.36 4.36 1 1];
        E=E*diag([rob,rob,rob,4.36*[1 1 1 1 1 1 1 1]]);
        rob=[1 1 1 1 1 1 4.36 4.36];
        X=X*diag([rob,rob,rob,[1 1 1 1 1 1 1 1]]);
        X=preproa(X+E,1);
    case 6,
        E=preproa(randn(n,32),1);
        %rob=[0 0 0 0 0 0 1 1];
        %X=X*diag([rob,rob,rob,[0 0 0 0 0 0 0 0]]);
        X=preproa(X+E,1);
    case 7,
        Xmean=inv(U'*U)*U'*X;
        B=orth(rand(8,2));
        C=orth(rand(4,2));
        [TB,TC,G]=tuck2a(Xmean,B,C,1);
        X=10*(X-U*Xmean)+U*G*kron(TC,TB)';
        %X=U*G*kron(TC,TB)'+2*randn(n,nv);    
    end
%
%X=X*diag(1./sqrt(diag(X'*X)/n));
%100*sum((U*pinv(U'*U)*U'*X-X).^2)./sum(X.^2)
%X(:,round(2*nv/3):nv)=X(:,round(2*nv/3):nv)*1.6;