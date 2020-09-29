function [f]=FunSEM(x)
global So nEXMV nENMV nEXLV nENLV pA pC pB pG pSF pSZ pSD pSE AA CC BB GG SF SD SE SZ Sp f Method

nMV=nEXMV+nENMV;
nLV=nEXLV+nENLV;
ntV=nEXLV+nENLV+nMV;

AA=zeros(nEXMV,nEXLV);
CC=zeros(nENMV,nENLV);
GG=zeros(nENLV,nEXLV);
BB=zeros(nENLV);
SF=zeros(nEXLV);
SZ=zeros(nENLV);
SD=zeros(nEXMV);
SE=zeros(nENMV);
m=1;
if pG(1,1)~=0
    for i=1:size(pG,1)
        GG(pG(i,1),pG(i,2))=x(m);
        m=m+1;
    end
end

if pB(1,1)~=0
    for i=1:size(pB,1)
        BB(pB(i,1),pB(i,2))=x(m);
        m=m+1;
    end
end

if pA(1,1)~=0
    for i=1:size(pA,1)
        AA(pA(i,1),pA(i,2))=x(m);
        m=m+1;
    end
    if size(AA,2)==1
        AA(find(AA==0),1)=1;
    else
        AA(find(sum(AA')==0),1)=1;
    end
end
if pC(1,1)~=0
    for i=1:size(pC,1)
        CC(pC(i,1),pC(i,2))=x(m);
        m=m+1;
    end
    if size(CC,2)==1
        CC(find(CC==0),1)=1;
    else
        CC(find(sum(CC')==0),1)=1;
    end
end
SF=eye(nEXLV);
if pSF(1,1)~=0
    for i=1:size(pSF,1)
        SF(pSF(i,1),pSF(i,2))=x(m);
        SF(pSF(i,2),pSF(i,1))=x(m);
        m=m+1;
    end
end
SZ=eye(nENLV);
if pSZ(1,1)~=0
    for i=1:size(pSZ,1)
        SZ(pSZ(i,1),pSZ(i,2))=x(m);
        SZ(pSZ(i,2),pSZ(i,1))=x(m);
        m=m+1;
    end
end

SD=eye(nEXMV);
if pSD(1,1)~=0
    for i=1:size(pSD,1)
        SD(pSD(i,1),pSD(i,2))=x(m);
        SD(pSD(i,2),pSD(i,1))=x(m);
        m=m+1;
    end
end

SE=eye(nENMV);
if pSE(1,1)~=0
    for i=1:size(pSE,1)
        SE(pSE(i,1),pSE(i,2))=x(m);
        SE(pSE(i,2),pSE(i,1))=x(m);
        m=m+1;
    end
end



F=[zeros(nEXLV,ntV);
GG BB zeros(nENLV,ntV-nEXLV-nENLV);
AA zeros(nEXMV,ntV-nEXLV);
zeros(nENMV,nEXLV) CC zeros(nENMV,ntV-nENLV-nEXLV)];


P=[SF zeros(nEXLV,ntV-nEXLV);
zeros(nENLV,nEXLV) SZ zeros(nENLV,ntV-nEXLV-nENLV);
zeros(nEXMV,nEXLV+nENLV) SD zeros(nEXMV,ntV-nEXMV-nENLV-nEXLV);
zeros(nENMV,ntV-nENMV) SE];

GF=[zeros(nMV,nLV) eye(nMV)];                   %filtering matrix
df = sum(1:size(GF,1)) - size(x,1);             % degrees of freedom
I=eye(ntV);	% Identidy matrix
Imv=eye(nMV);
Sp = GF * inv(I-F) * P * (inv(I-F))' * GF';     % PREDICTED covariance (or correlation) matrix

						% alternative fitting functions
switch Method
   case 'ML'
      f = log(det(Sp))-log(det(So)) + trace(So*(pinv(Sp)))-nMV;
   case 'OLS'
      f = (trace((So-Sp)*(So-Sp)));
   case 'GLS'
      f = 0.5 * trace((Imv-(pinv(Sp)*So))*(Imv-(pinv(Sp)*So)));
    case 'WLS'
      f = 0.5 * trace((Imv - (pinv(So)*Sp))*(Imv - (pinv(So)*Sp)));
   otherwise
      disp('Error: Unknown method')
end
                     



