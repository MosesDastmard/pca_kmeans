%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEM Stima sequenziale         %
% soluzione dei minimi quadrati %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IntelligentData
%Y=X(:,5:7);
%X=X(:,1:4);
n=100;
J=eye(n)-(1/n)*ones(n);
Sx=(1/n)*X'*J*X;
Sy=(1/n)*Y'*J*Y;
Xs=J*X*diag(diag(Sx))^-0.5;
Ys=J*Y*diag(diag(Sy))^-0.5;
[A,psi,T,stats,F] = factoran(Xs,1);
[C,psi2,T2,stats2,N] = factoran(Ys,1);
%[At,Ft,eige1] = princomp(Xs);
%A=At(:,1);
%F=Xs*A;
%[Ct,Nt,eige2] = princomp(Ys);
%C=Ct(:,1);
%N=Ys*C;
SF=(1/n)*F'*F;
SN=(1/n)*N'*N;
D=X-F*A';
SD=(1/n)*D'*D;
E=Y-N*C';
SE=(1/n)*E'*E;
G=pinv(F)*N;
B=0;
Z=N-F*G;
SZ=(1/n)*Z'*Z;
SSxy=A*SF*G'*C';
SSx=A*SF*A'+SD;
SSy=C*SN*C'+SE;
So=(1/n)*[X Y]'*J*[X Y];
SS=[SSx SSxy; SSxy' SSy];
f= trace((So-SS)'*(So-SS))
%  step2

N=F*G;
SN=(1/n)*N'*N;
C=(pinv(N)*Y)';
E=Y-N*C';
SE=(1/n)*E'*E;
G=pinv(F)*N;
Z=N-F*G;
SZ=(1/n)*Z'*Z;
SSxy=A*SF*G'*C';
%SSx=A*SF*A'+SD;
SSy=C*SN*C'+SE;
SS=[SSx SSxy; SSxy' SSy];
f= trace((So-SS)'*(So-SS))


