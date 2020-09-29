function [Tab,GLM,GLM_R2,GLM_AIC,GLM_BIC,ARI,MIS,PCHI2]=evasim(n, K, Pr, nsim, esigma)

%Initialization indices
ARI = zeros(nsim,1);
Tab = zeros(K,K,nsim);
GLM_R2 = zeros(nsim,1);
GLM_AIC = zeros(nsim,1);
GLM_BIC = zeros(nsim,1);
PCHI2 = zeros(nsim,1);
MIS = zeros(nsim,1);

for sim=1:nsim
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a mixture of Gaussian distributions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist=6;
csigma=0;

[X,U,u]=MixSampling(n, Pr, dist, csigma, esigma);


[loopOtt,UOtt,fOtt,iterOtt,cluster]=kmeansVICHI(X,K,10); %% Application K-Means

% Ordinal regression model 

%GLM = fitglm(X,cluster,'distr','normal');                                 %logistic regression
GLM = stepwiseglm(X,cluster,'constant','upper','linear','distr','normal'); %stepwise  
fit = round(GLM.Fitted.Response);

fitted = zeros(n,K);

for i=1:n
    for k=1:K
        if fit(i,1)==k
            fitted(i,k)=1;
        end
    end
end

[tbl,chi2,p] = crosstab(cluster, fit);

% solve lable switching problem

C=zeros(K);
    for k=1:K
        for g=1:K
             C(k,g)=sum((UOtt(:,k)-fitted(:,g)).^2 );
        end
    end
[ss,unfeas,DD] = bghungar(-C);
fitted = fitted(:,ss);
table = UOtt'*fitted;

% saving
Tab(:,:,sim) = table;
ARI(sim,1) = mrand(table);
MIS(sim,1) = 1-sum(diag(table))/sum(sum(table));
GLM_R2(sim,1) = GLM.Rsquared.Adjusted;
GLM_AIC(sim,1) = GLM.ModelCriterion.AIC;
GLM_BIC(sim,1) = GLM.ModelCriterion.BIC;
PCHI2(sim,1) = p;

   if sim <= nsim
        disp(sprintf('n.Random=%g, GLM_R2=%g, GLM_AIC=%g, GLM_BIC=%g, ARI=%g, MIS=%g, PCHI2=%g',sim, GLM_R2(sim), GLM_AIC(sim), GLM_BIC(sim), ARI(sim), MIS(sim), PCHI2(sim)))
   end
end

%plot(X(1:100,1), X(1:100,2), 'or')
%hold on
%plot(X(101:200,1), X(101:200,2), '+b')
%plot(X(201:300,1), X(201:300,2), '*g')

%plotmatrix(X);

figure
subplot(3,1,1)
hist(ARI,40)
title('Adjusted Rand Index')
xlim([0 1])
subplot(3,1,2)
hist(GLM_R2,40)
title('R squared index by Generalized Linear Model')
xlim([0 1])
subplot(3,1,3)
hist(MIS,40)
title('MisClassification rate')
xlim([0 1])


% Related functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-means algorithm di VICHI%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [loopOtt,UOtt,fOtt,iterOtt,cluster]=kmeansVICHI(X,K,Rndstart)
%
% n = number of objects
% J = number of variables
% K = number of clusters of the partition
%
% maxiter=max number of iterations
%

maxiter=100;
n = size(X,1);
J = size(X,2);
epsilon=0.000001;


% initial partition U0 is given

% best in a fixed number of partitions
%seed=200;             %si può rimettere per ritrovare le soluzioni ottenute
%rand('state',seed)    %
for loop=1:Rndstart
   U0=randPU(n,K);
%
   su=sum(U0);
   %
   % given U compute Xmean (compute centroids)
   Xmean0 = diag(1./su)*U0'*X;

   for iter=1:maxiter
   %
   % given Xmean0 assign each units to the closest cluster
   %
        U=zeros(n,K);
        for i=1:n
            mindif=sum((X(i,:)-Xmean0(1,:)).^2);
            posmin=1;
            for j=2:K
                dif=sum((X(i,:)-Xmean0(j,:)).^2);
                if dif < mindif
                    mindif=dif;
                    posmin=j;
                end 
            end
            U(i,posmin)=1;
        end
   % given a partition of units 
   % i.e, given U compute Xmean (compute centroids)
   %
        su=sum(U);
        while sum(su==0)>0,
            [m,p1]=min(su);
            [m,p2]=max(su);
            ind=find(U(:,p2));
            ind=ind(1:floor(su(p2)/2));
            U(ind,p1)=1;
            U(ind,p2)=0;
            su=sum(U);
        end 
   %
   % given U compute Xmean (compute centroids)
        Xmean = diag(1./su)*U'*X;
   
   %
   % 
   % compute objective function
   %
        BB=U*Xmean-X;
        f=trace(BB'*BB);
%   
%  stopping rule
%   
        dif=sum(sum((Xmean-Xmean0).^2));
        if dif > epsilon  
            Xmean0=Xmean;
        else
            break
        end
   end
   if loop==1
        UOtt=U;
        fOtt=f;
        loopOtt=1;
        iterOtt=1;
        XmeanOtt=Xmean;
   end
   if f < fOtt
        %disp(sprintf('k-means: loopOtt=%g, fOtt=%g, iter=%g',loopOtt,fOtt,iter))
        UOtt=U;
        fOtt=f;
        loopOtt=loop;
        iterOtt=iter;
        XmeanOtt=Xmean;
   end
   
end
        %disp(sprintf('k-means: loopOtt=%g, fOtt=%g, iter=%g',loopOtt,fOtt,iter))

   cluster = zeros(n,1);
   
   for t=1:n
       for s=1:K
           if U(t,s)==1
           cluster(t,1)=s;
           end
       end
   end

function [U]=randPU(n,c)

% generates a random partition of n objects in c classes
%
% n = number of objects
% c = number of classes
%
U=zeros(n,c);
U(1:c,:)=eye(c);

U(c+1:n,1)=1;
for i=c+1:n
    U(i,[1:c])=U(i,randperm(c));
end
U(:,:)=U(randperm(n),:);



function mri=mrand(N)
%
% modified rand index (Hubert & Arabie 1985, JCGS p.198)
%
n=sum(sum(N));
sumi=.5*(sum(sum(N').^2)-n);
sumj=.5*(sum(sum(N).^2)-n);
%pb=sumi*sumj/(n*(n-1));
%mri=(.5*(sum(sum(N.^2))-n)-pb)/((sumi+sumj)/2-pb);

% nuova versione
pb=sumi*sumj/(n*(n-1)/2);
mri=(.5*(sum(sum(N.^2))-n)-pb)/((sumi+sumj)/2-pb);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a mixture of Gaussian distributions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,U,u]=MixSampling(n, Pr, dist, csigma, esigma)
% input
% n: numebr of Objects; 
% Pr: vector of a priori probabilities (K x 1);
% dist: constant distance in symmetric matrix; 
% csifma: error to perturbate centroids;
% esigma: standard deviation of error

%output
% X Data matrix (n x K-1)
% U membership matrix
% u categorical clustering variable 

% K number of cluster
K=size(Pr,1);

% Compute Xm: K equidistant centroid matrix
D=dist*(ones(K)-eye(K));
Jc=eye(K)-(1/K)*ones(K);
cDc=-0.5*Jc*D*Jc;
[A,B]=eigs(cDc,K-1);
Xm=A*B.^0.5;
%
%
% Compute centroid matrix
%
mu=zeros(K-1,1);
Sigma=1*eye(K-1);
E = mvnrnd(mu,Sigma,K)*csigma;
%
% add an error to centroids (modify isolation)
%
Xmp=Xm+E;
%
% variance and covariance 
%
SigmaO = cov(rand(100,K-1),1)*esigma; % modify heterogeneity

% Compute n multinomial observations (1 unit sample)

U = multrnd(1,Pr,n);
u=U*[1:K]';
[u,inu] = sort(u);
U=U(inu,:);
for i=1:n
    X(i,:)=mvnrnd(Xmp(u(i),:),SigmaO,1);
end

if K<=8

axis square
gplotmatrix(X,[],u,'brkmcywv','+x^d.*os',[],'off','hist');

end


% ************************************************
%   see my contrib. to file exchange: this is an improved (more robust)
%   version.

function [iz, unfeas, D] = bghungar(C)
%BGHUNGAR "Hungarian algorithm" to solve the square assignment problem
%
%#       For:
%#        C - a square profit/cost matrix.
%  [iz, D] = bghungar(C);
%#       Returns:
%#       iz - the optimal assignment: MAXIMIZES total profit
%#        D - the square matrix equivalent to C at the end of iteration [1]
%
% See also the script T.M containing an example matching N=6 pairs of signals
% from 2 experiments. In it the A(!)symmetric matrix C contains the cross
% correlation coefficients for each possible pairwise matching
%  Note: For assignments that MINIMIZE cost just call BGHUNGAR, inverting
%        the sign of the cost matrix
%         [iz, D] = bghungar(-C);

%# (c) 2002/11/29 Nedialko Krouchev, Krouchen@physio.umontreal.ca
%# Comments and feedback on code quality issues very welcome
%# original & pure Matlab implementation
%# Copyright (c) 1989-2002 Nedialko Krouchev
%# $Revision: 1.1 $  $Date: 29-Nov-2002 11:43:33 $

%# This code is made available according to conditions in the general spirit
%# of the FSF/GNU public license; it is freeware, but not shareware - the author
%# reserves all rights over its intellectual content
%# Applicability and proper use of this code is entirely each user's responsibility.

%# References:
%# [1] D.Ivantchev and G.Negler "Network Optimization", Sofia 1992

%# ===========================================================================

unfeas = 0;

%#(1)
[m,n] = size(C);

[r,ii] = max(C,[],1); C1 = ones(m,1)*r - C;

[c,jj] = min(C1,[],2); D = C1 - c*ones(1,n);

iz = zeros(m,1);
if any(isnan(D)), unfeas = 1; return; end

for j = 1:m
   kk = find( ~D(:,j) );
   while ~isempty(kk),
      i=kk(1); if ~iz(i), iz(i) = j; break; else, kk(1) = []; end
   end
end %j

izPrev = zeros(1,m);
%# --------------------------------------------

%#(2)
keepSets = 0;

%# occupied rows:
ii = find( iz );
while length( ii ) < m


if ~keepSets
   keepSets = 0;

%# ============================================
%# Reset the free sets:

   %# all rows & columns:
   rr = 1:m; cc = 1:n;

   %# --------------------------------------------

   %#(3)

   %# occupied columns:
   zz = iz(ii);

   %# free rows & columns:
   %% rr( ii ) = [];
   cc( zz ) = [];

%# ============================================

end %if ~keepSets

%# --------------------------------------------

while 1
   if any(isnan(D)), unfeas = 2; break; end

%#(4,5)
%# Find a row containing free zeros:

jz = 0;
for i = rr
   kk = find( ~D(i,cc) );
   if ~isempty(kk), jz = i; kz = cc(kk(1)); break; end
end %i

%# --------------------------------------------

if jz
%# Found a row containing free zeros:

   if iz(jz)
%#(6)
%# This row already has an assigned zero:
      kz = iz(jz); cc = [cc,kz]; zz(zz==kz) = []; rr(rr==jz) = [];

   else

%#(7)
%# This row has no assigned zero yet:
%#"Chain": Find a way to assign a zero to it:

% if size( izPrev, 1 ) == 3, disp('Start debugging: dbstep'); keyboard; end
      jz0 = jz; kz0 = kz; iz0 = iz;

      while 1

%# in a column:
         iz(jz) = kz;
         rr1 = [1:jz-1, jz+1:m]; next = find( iz(rr1) == kz );
         if isempty(next), break; end
         jz = rr1(next(1));

%# in a row:
         iz(jz) = 0;
         cc1 = [1:kz-1, kz+1:n]; next = find( ~D(jz,cc1) );
         if isempty(next), break; end
         kz = cc1(next(1));

      end
%# go out of the infinite while-loop %#(5) here.
%# Detect alternating chains:
      if any( ismember( iz', izPrev,  'rows') )
         unfeas = 3; % break;

% kz0, iz', izPrev, iz0', disp(': dbcont'); keyboard
% jz0, kz0, rr, cc, E=D; E(logical(D))=1, disp(': dbstep'); keyboard

%# Block the free row/col (jz0,kz0) where all this started;
%# Recover last iz0 (iz0)  & reset izPrev ???

         iz = iz0;
         zz = [zz;kz0]; cc(cc==kz0) = []; %% rr(rr==jz0) = [];
% rr, cc

      else
        izPrev = [ izPrev; iz' ];
      end
      break;
   end

%# --------------------------------------------
else
%# No row containing free zeros found:
%# Create free zeros:

%#(8)
   p = min( min( D(rr,cc) ) );
   if p >= Inf, unfeas = 4; break; end
   D(rr,:) = D(rr,:) - p;
   D(:,zz) = D(:,zz) + p;

end% if jz

%# --------------------------------------------
end %while %#(5)

if unfeas==3,
   unfeas = 0; keepSets = 1;
else
   keepSets = 0;
end

if unfeas, break; end

ii = find( iz );

%# ============================================
end %while %#(2)


%# ============================================
% C, iz', izPrev, E=D; E(logical(D))=1, disp('Final: dbstep'); keyboard
% keepSets = 0;
