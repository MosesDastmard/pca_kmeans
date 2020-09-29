% Matrix B for CluReg_LS
% versione semplificata di PesiB
% 21 settembre 2007

function  B=PesiB2(X,Y,V,C)

[I,J]=size(X);
[M]=size(Y,2);
C_tilde=V*C;
MAT=zeros(M*I,J);
m=0;  
[Q,gf]=size(C);


HH=corr(X);
for q=1:Q
a=find(V(:,q));
[dim,se]=size(a);
HH([a],[a])=zeros(dim,dim);
end;
% 

   for j=1:J
       m=m+1;
       x_tildej=X(:,j); c_tildej=C_tilde(j,:)';  
       MAT(:,m)=kron(c_tildej,x_tildej);
   end
   b=pinv(MAT'*MAT)*MAT'*Y(:); 
   %b=pinv(HH+1*eye(J,J))*Y';
   B=diag(b);
%%%%%%%%%%%%%%%%%%%%%%%%%%


  