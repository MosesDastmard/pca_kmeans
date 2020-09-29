% generate data from Dppca model

 % X = YV'B + E 
 
 function [Xs,VG,BG,YG]=GenDataDPPCA(n,J,K,PercErr)
 A=zeros(J,K);
 VG=randPU(J,K);
 v=VG*[1:K]';
 [in,ln]=sort(v);
 VG=VG(ln,:);
 Y1=rand(n,K);
 Y1=zscore(Y1);
 Y=orth(Y1);
 D=diag(diag(cov(Y,1)).^0.5);
 YG=Y*D^-1;
 b=rand(1,J);
 for k=1:K
    iCk=find(VG(:,k)==1);
    sbk=sum(b(iCk));
    b(iCk)=sqrt(b(iCk)./sbk);
 end
 BG=diag(b);
 E=randn(n,J);
 X=YG*VG'*BG+PercErr*E;
 Xs=zscore(X,1);
 
 s2=0;
 for g=1:K
    ibCg=VG(:,g); 
    JCg=find(ibCg==1);
    Sg=(1./n)*Xs(:,JCg)'*Xs(:,JCg);
    [U,L]=eig(Sg);
    [l,il]=sort(diag(L), 'descend');
    L=diag(l);
    U=U(:,il);
    A(JCg,g)=U(:,1)*l(1).^0.5;
    [l,il]=sort(diag(L), 'descend');
    s2 = s2+sum(l(2:size(Sg,1)));                     
 end
 s2=s2./(J-K);
 %A=BG*VG;
 S=(1/n)*Xs'*Xs;
 Sx = A*A' + s2*eye(J);
 f = (-n/2)*( J*log(pi) + log(det(Sx)) + trace(Sx\S) )

 
 
 
 
 