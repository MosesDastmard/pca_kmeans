%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structural equation model SEM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maurizio Vichi
% November 2012

function [AA,CC,GG,BB, SF, SZ,SD, SE,  f]=SEM(X, Y,td, mpA,mpC,mpB,mpG,mpSF,mpSZ,mpSD,mpSE, Mt)
global So nEXMV nENMV nEXLV nENLV pA pC pB pG pSF pSZ pSD pSE AA CC BB GG SF SD SE SZ Sp f Method
% method of estimation
Method=Mt;
pA=mpA;
pB=mpB;
pC=mpC;
pG=mpG;
pSF=mpSF;
pSZ=mpSZ;
pSD=mpSD;
pSE=mpSE;
[n,J]=size(X);
[M]=size(Y,2);
% parameters to estimate

npA=size(pA,1);
if npA==1 & pA(1,1)==0
    npA=0;
    nEXLV=0;
else
    [nEXLV]=max(pA(:,2));
end
npC=size(pC,1);
if npC==1 & pC(1,1)==0
    npC=0;
    nENLV=0;
else
    [nENLV]=max(pC(:,2));
end
if nENLV==0
    nENLV=max(pG(:,1));
end
npB=size(pB,1);
if npB==1 & pB(1,1)==0
    npB=0;
else
    npB=size(pB,1);
end
npG=size(pG,1);
if npG==1 & pG(1,1)==0
    npG=0;
else
    npG=size(pG,1);
end
if pSF(1,1)==0
    npSF=0;
else
    npSF=size(pSF,1);
end
if pSZ(1,1)==0
    npSZ=0;
else
    npSZ=size(pSZ,1);
end
if pSD(1,1)==0
    npSD=0;
else
    npSD=size(pSD,1);
end
if pSE(1,1)==0
    npSE=0;
else
    npSE=size(pSE,1);
end
nEXMV=J;
nENMV=M;
nMV=nEXMV+nENMV;
nLV=nEXLV+nENLV;
ntV=nEXLV+nENLV+nMV;
np=npA+npC+npB+npG+npSF+npSZ+npSD+npSE;
XY=[X,Y];
Jm=eye(n)-(1/n)*ones(n);
% compute covariance matrix of X and Y
Sx=(1/n)*X'*Jm*X;
if M~=0
    Sy=(1/n)*Y'*Jm*Y;
end
if td == 'S'
% Standardize X and Y
    Xs=Jm*X*diag(diag(Sx).^-0.5);
    if M~=0
        Ys=Jm*Y*diag(diag(Sy).^-0.5);
    else
        Ys=[];
    end
    XYs=[Xs Ys];
% compute covariance of standardized data [Xs Ys] 
    So=(1/n)*(XYs'*XYs);
else
% compute covariance of non-standardized data[X Y]    
    So=(1/n)*(XY'*Jm*XY);
end
% Starting  solution Sequential Estimation of SEM

%
%[a,psi,T,stats,F] = factoran(Xs,1);
%[c,psi2,T2,stats2,N] = factoran(Ys,1);

A1=[];
b1=[];
A0=[];
b0=[];
nonlcon=[];
if td=='S'
    lb=-1*ones(nMV+npB+npG,1);
    lb(nMV+npB+npG+1:np)=0.001*ones(np-nMV-npB-npG,1);
    lu=1*ones(np,1);
else
    lb=-10000000000000*ones(nMV+npB+npG,1);
    lb(nMV+npB+npG+1:np)=0.001*ones(np-nMV-npB-npG,1);
    lu=10000000000000*ones(np,1);
end
options=optimset('Algorithm','sqp','Display','off');
x=2*(0.5-rand(nMV+npB+npG,1));
x(nMV+npB+npG+1:np)=ones(np-nMV-npB-npG,1); 

x=fmincon(@FunSEM,x,A1,b1,A0,b0,lb,lu, nonlcon, options);



						% alternative fitting functions
I=eye(nMV);					                     % Identidy matrix
df=nMV*(nMV+1)/2-np;    % degree of freedom
X2=(n-1)*f;
pValue =1 - chi2cdf(X2, df);
piSp=pinv(Sp);
switch Method
   case 'ML'
      disp('SEM ESTIMATION: Maximum Likelihood')
      % compute Chi-square test
      GFI= 1 - trace((piSp*So - I)*(piSp*So - I))/trace((piSp*So)*(piSp*So));
    case 'OLS'
      disp(sprintf('\n SEM ESTIMATION: Ordinary Least Squares \n' ))
      GFI=1-trace((So-Sp)*(So-Sp))/trace(So*So);
   case 'GLS'
      disp(sprintf('\n SEM ESTIMATION: Generalized Least Squares \n'))
      GFI= 1- trace((I-(piSp*So))*(I-(piSp*So)))/nMV;
    case 'WLS'
      disp(sprintf('\n SEM ESTIMATION: Weighted Least Squares \n'))
      GFI = trace((I - (inv(So)*Sp))*(I - (inv(So)*Sp)))/nMV;
   otherwise
      disp('\n Error: Unknown method')
end
AGFI=1 - (nMV*(nMV+1)/(2*df))*(1-GFI);
RMSEA = sqrt(max([f /(df)-1/(n-1), 0]));

disp(sprintf('Number of Cases  = %g', n))
disp(sprintf('Number of Exogenous  Latent   Variables = %g ', nEXLV))
disp(sprintf('Number of Endogenous Latent   Variables = %g', nENLV))
disp(sprintf('Number of Exogenous  Manifest Variables = %g', nEXMV))
disp(sprintf('Number of Endogenous Manifest Variables = %g \n', nENMV))

switch td
    case 'S'
        disp('Scale of Data = Standardized Data')
    case 'N'
        disp('Scale of Data = No Standardized Data \n')
end 

disp(sprintf('Number of Observed Statistics = %g', nMV*(nMV+1)/2))
disp(sprintf('Number of Estimated Parameters = %g \n \n', np))


disp(sprintf('Model Chi-Square Test Chi2-value =%g  Degree of freedom=%g p-value=%g \n', X2,df,pValue ))
disp(sprintf('FIT of the STRUCTURAL EQUATION MODEL \n Goodness-of-FIT  GFI=%g \n Adjusted GFI   AGFI=%g', GFI,AGFI))
disp(sprintf('\n Root-Mean-Square-Error-of-Approximation \n RMSEA=%g', RMSEA))


disp(sprintf('\n \n Measurement Model for Exogenous Variables \n'))

disp(sprintf('\n Parameter Estimates \n \n \n'))
disp(sprintf('     Path          Path coef       Std Error        Pr(p>|Z|)        Var Error        Communality' ))

for h=1:nEXLV
    p=zeros(J,1);
    %p=zeros(J,2);
    for j=1:size(pA,1)
        %p(j,:)= normcdf([-1*abs(AA(j,1)/SD(j,j).^0.5) abs(AA(j,1)/SD(j,j).^0.5)]);
        if AA(j,h)~=0
            p(j)= r2pv(AA(j,h),n);
            dSD=diag(SD)./n;
            disp(sprintf('X(%g) <-- Xi(%g)     % f        %f         %f         %f         %f', pA(j,1), pA(j,2), AA(j,h), (dSD(j)).^0.5, p(j), dSD(j).*n, AA(j,h).^2))
    
        end
    end
end
if nEXLV>1
    disp(sprintf('\n \n Correlation Matrix of Latent Exogenous Variables'))
    SF
end



if size(CC,1)~=0 
    disp(sprintf('\n \n Measurement Model for Endogenous Variables \n'))
    disp(sprintf('     Path           Path coef       Std Error        Pr(p>|Z|)        Var Error        Communality' ))

    for h=1:nENLV
        p=zeros(M,1);
        for j=1:size(pC,1)
            if CC(j,h)~=0
                p(j)= r2pv(CC(j,h),n);   
                dSE=diag(SE)./n;
                disp(sprintf('Y(%g) <-- Eta(%g)     % f        %f         %f         %f         %f', pC(j,1), pC(j,2), CC(j,h), (dSE(j)).^0.5, p(j), dSE(j).*n, CC(j,h).^2))
            end
        end
    end
end
if sum(sum(GG))~=0
    disp(sprintf('\n \n Structural Model Exogenous Variables\n'))
    disp(sprintf('\n Parameter Estimates \n  \n \n'))
    disp(sprintf('     Path            Path coef       Std Error        Pr(p>|Z|)' ))
    for h=1:size(pG,1)
        pp = r2pv(  GG(pG(h,1),pG(h,2)),n);   
        dSD=diag(SD)./n;
        disp(sprintf('Eta(%g) <-- Xi(%g)     % f        %f         %f         %f', pG(h,1),pG(h,2),GG(pG(h,1),pG(h,2)), dSD(pG(h,1)).^0.5, pp))
     end
end


if sum(sum(BB))~=0
    disp(sprintf('\n \n Structural Model Endogenous Variables\n'))
    disp(sprintf('\n Parameter Estimates \n  \n \n'))
    disp(sprintf('     Path             Path coef        Std Error        Pr(p>|Z|)        Var Error' ))

    for h=1:size(pB,1)
        pp = r2pv(BB(pB(h,1),pB(h,2)),n);   
        dSZ=diag(SZ)./n;
        disp(sprintf('Eta(%g) <-- Eta(%g)     % f        %f         %f         %f', pB(h,1), pB(h,2),BB(pB(h,1),pB(h,2)), dSZ(pB(h,1)).^0.5, pp, dSZ(pB(h,1)).*n))
     end
end
    
    




 function p=r2pv(r,n)
%
% 	p=r2pv(r,n)
%
% r = estimated correlation coefficient (IE |r| <= 1)
%   = (1/n)*(x'*y) for col vectors x,y of length n
% n = no. samples used
% p = P-value based on |r| (two sided) with rho=0 (null case)
%
% NOTES: following Cramer, p.400, convert r to a t and use what we have for t 
if n < 3
    error('n < 3');
end
if r==1. 
    p=0; 
    return;
end
t=sqrt(n-2)*r/(sqrt(1-r*r)); 	% this is t with n-2 d.f.
t=abs(t);							% use |t| for two sided P-value
p=2*(1-tcdf(t,n-2));





