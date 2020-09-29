%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§%%%%%%%
% Algorithm for Maximum Likelihood Mixture Clustering  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% for PhD 2008
%
function [loopOtt,UOtt,LLikeOtt,pOtt,iterOtt]=MLMC(X,K,Rndstart)
%
% I = number of objects
% J = number of variables
% K = number of clusters of the partition
%
% maxiter=max number of iterations
%
% Rndstart = number of initial random starts
%
maxiter=100; % 
I = size(X,1);
J = size(X,2);
epsilon=0.000000000001;

%seed=2;             %you can restart with a fixed seed of the random numbers
%rand('state',seed)    %
%
for loop=1:Rndstart
    
   % initialize parameters
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Initial partition U0 is given  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U0=randPU(I,K);
   %
   su=sum(U0);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given U, compute initial mean vector (compute centroids)   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Xmean0 = diag(1./su)*U0'*X;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given mean vector and U0, compute initial Var-Covar Matrix %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Sig0=zeros(J);
   for k=1:K
      for i=1:I,
         Sig0=Sig0+U0(i,k)*((X(i,:)-Xmean0(k,:))'*(X(i,:)-Xmean0(k,:)));
      end
   end
   Sig0=(1/I)*Sig0;
   Sig0m1=inv(Sig0);
   detSig0=det(Sig0);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   % given U0, compoute the initial prior probabilities  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   p0=(1/I)*sum(U0)';
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % compute intitial Log-Likelihhod
   %
   LLike0=0;
   for i=1:I
      for k=1:K
                LLike0=LLike0 + U0(i,k)*((-(J./2)*log(2*pi))-(0.5*log(detSig0))-0.5*((X(i,:)-Xmean0(k,:))*Sig0m1*(X(i,:)-Xmean0(k,:))')+ log(p0(k)));
      end
   end
   entra=0; 
   for iter=1:maxiter
 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % given Xmean0, assign each units to the closest Gaussian component %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
        Sig0m1=inv(Sig0);
        detSig0=det(Sig0);
        U=zeros(I,K);
        for i=1:I
            %mindif =(X(i,:)-Xmean0(1,:))*Sig0m1*(X(i,:)-Xmean0(1,:))' - log(p0(1));
            mindif=(((J./2)*log(2*pi))+(0.5*log(detSig0))+0.5*((X(i,:)-Xmean0(1,:))*Sig0m1*(X(i,:)-Xmean0(1,:))')- log(p0(1)));
            posmin=1;
            for k=2:K
                %dif=(X(i,:)-Xmean0(k,:))*Sig0m1*(X(i,:)-Xmean0(k,:))' - log(p0(k));
                dif=(((J./2)*log(2*pi))+(0.5*log(detSig0))+0.5*((X(i,:)-Xmean0(k,:))*Sig0m1*(X(i,:)-Xmean0(k,:))')- log(p0(k)));
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
            entra=1;        
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
        % given U, compute Xmean (compute centroids) and Variance-Covariance matrix %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        su=sum(U);        
        Xmean = diag(1./su)*U'*X;
        
        Sig=zeros(J);
        for k=1:K
            for i=1:I,
                Sig = Sig + U(i,k)*((X(i,:)-Xmean(k,:))'*(X(i,:)-Xmean(k,:)));
            end
        end
        Sig=(1./I)*Sig;
        Sigm1=inv(Sig);
        detSig=det(Sig);

        %      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % given U, compute the prior probabilities   %   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        p=(1./I)*sum(U)';

%
% 
% compute Log-Likelihood function
   %
        LLike=0;
        for i=1:I
            for k=1:K
                LLike = LLike + U(i,k)*((-(J./2)*log(2*pi))-(0.5*log(detSig))-0.5*((X(i,:)-Xmean(k,:))*Sigm1*(X(i,:)-Xmean(k,:))')+ log(p(k)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%   
        %       stopping rule     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %
        dif=LLike-LLike0;
       %
       %
        if dif < 0 & entra==0
            disp(sprintf('MaxLikeMixClus: Log-Likelihood decreases iter=%g, LLike0=%g, dif=%g', iter, LLike0, dif))
        end
        if dif > epsilon  
            Xmean0=Xmean;
            Sig0=Sig;
            LLike0=LLike;
            p0=p;
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
        pOtt=p;
   end
   if dif < 0 & entra==1
   disp(sprintf('MaxLikeMixClus: loop=%g, LLike=%g, iter=%g warning convergence',loop,LLike,iter))
   else
   disp(sprintf('MaxLikeMixClus: loop=%g, LLike=%g, iter=%g',loop,LLike,iter))
   end

   if LLike > LLikeOtt
        UOtt=U;
        LLikeOtt=LLike;
        loopOtt=loop;
        iterOtt=iter;
        XmeanOtt=Xmean;
        pOtt=p;
   end
end
disp(sprintf('MaxLikeMixClus: loopOtt=%g, fOtt=%g, iter=%g',loopOtt,LLikeOtt,iter))