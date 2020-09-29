%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm for Maximum Likelihood Clustering  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% for PhD 2008
%
function [loopOtt,UOtt,LLikeOtt,iterOtt]=MLC(X,K,Rndstart)
%
% I = number of objects
% J = number of variables
% K = number of clusters of the partition
%
% maxiter=max number of iterations
%
%
% Rndstart = number of initial random starts
%

maxiter=100; % 
I = size(X,1);
J = size(X,2);
epsilon=0.000001;

% initial partition U0 is given

% best in a fixed number of partitions

%seed=1;             %you can restart with a fixed seed of the random numbers
%rand('state',seed)    %
for loop=1:Rndstart
    %
    % Initialize parameters
    %
    U0=randPU(I,K);
%
    su=sum(U0);
   %
   % given initial U, compute initial mean vector (compute centroids)
   %
   Xmean0 = diag(1./su)*U0'*X;
   %
   % given intitial  mean vector and U0, compute initial Var-Covar Matrix
   %
   % compute intitial Log-Likelihhod
   %
   Sig0=zeros(J);
   for k=1:K
      for i=1:I,
         Sig0=Sig0+U0(i,k)*((X(i,:)-Xmean0(k,:))'*(X(i,:)-Xmean0(k,:)));
      end
   end
   Sig0=(1/I)*Sig0;
   Sig0m1=inv(Sig0);
   LLike0=0;
   for i=1:I
      for k=1:K
                LLike0=LLike0 + U0(i,k)*((-(J/2)*log(2*pi))-(0.5*log(det(Sig0)))-0.5*((X(i,:)-Xmean0(k,:))*Sig0m1*(X(i,:)-Xmean0(k,:))'));
      end
   end

   for iter=1:maxiter
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given Xmean0 assign each units to the closest cluster %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
        U=zeros(I,K);
        for i=1:I
            mindif =(X(i,:)-Xmean0(1,:))*inv(Sig0)*(X(i,:)-Xmean0(1,:))';
            posmin=1;
            for k=2:K
                dif=(X(i,:)-Xmean0(k,:))*inv(Sig0)*(X(i,:)-Xmean0(k,:))';;
                if dif < mindif
                    mindif=dif;
                    posmin=k;
                end 
            end
            U(i,posmin)=1;
        end
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
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given U, compute Xmean (compute centroids) ans Variance-Covariance matrix %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
        Xmean = diag(1./su)*U'*X;

        Sig=zeros(J);
        for k=1:K
            for i=1:I,
                Sig=Sig+U(i,k)*((X(i,:)-Xmean(k,:))'*(X(i,:)-Xmean(k,:)));;
            end
        end
        Sig=(1/I)*Sig;
        Sigm1=inv(Sig);

   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % compute Log-Likelihood function %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
        LLike=0;
        for i=1:I
            for k=1:K
                LLike=LLike + U(i,k)*((-(J/2)*log(2*pi))-(0.5*log(det(Sig)))-0.5*((X(i,:)-Xmean(k,:))*Sigm1*(X(i,:)-Xmean(k,:))'));
            end
        end
 
   %   
   %%%%%%%%%%%%%%%%%%%   
   %  stopping rule  %
   %%%%%%%%%%%%%%%%%%%
   %   
        dif=LLike-LLike0;
    if dif < 0
        disp(sprintf('MaxLikeClus: ATTENTION! - Log-Likelihood decreases'))
    end
        if dif > epsilon  
            Xmean0=Xmean;
            Sig0=Sig;
            LLike0=LLike;
        else
            break
        end
   end
   if loop==1
        UOtt=U;
        LLikeOtt=LLike;
        loopOtt=1;
        iterOtt=1;
        XmeanOtt=Xmean;
   end
   disp(sprintf('MaxLikeClus: loop=%g, LLike=%g, iter=%g',loop,LLike,iter))
   if LLike > LLikeOtt
        UOtt=U;
        LLikeOtt=LLike;
        loopOtt=loop;
        iterOtt=iter;
        XmeanOtt=Xmean;
   end
end
 disp(sprintf('MaxLikeClus: loopOtt=%g, LLikeOtt=%g, iterOtt=%g',loopOtt,LLikeOtt,iterOtt))


