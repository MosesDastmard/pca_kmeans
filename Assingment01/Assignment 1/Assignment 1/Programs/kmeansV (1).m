%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-means algorithm for data matrices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vichi Maurizio 
% version 10.10.2012

function [loopOtt,UOtt,fOtt,iterOtt]=kmeansV(X,K,Rndstart)
%
% n = number of objects
% J = number of variables
% K = number of clusters of the partition
%
% maxiter=max number of iterations
% X data Matrix
% Rndstart = number of random starts
% 
% initialization
%
maxiter=100;
n = size(X,1);
J = size(X,2);
epsilon=0.000001;

% find the best solution in a fixed number of random start partitions
%seed=200;            
%rand('state',seed)    %
for loop=1:Rndstart
    
   % initial partition U0 is given
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
        disp(sprintf('k-means: loopOtt=%g, fOtt=%g, iter=%g',loopOtt,fOtt,iter))
        UOtt=U;
        fOtt=f;
        loopOtt=loop;
        iterOtt=iter;
        XmeanOtt=Xmean;
   end
end
        disp(sprintf('k-means: loopOtt=%g, fOtt=%g, iter=%g',loopOtt,fOtt,iter))

