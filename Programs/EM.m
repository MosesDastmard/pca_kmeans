%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§%%%%%%%
% EM Algorithm Gaussian Mixture Clustering             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% for PhD 2008
%
function [loopOtt,UOtt,LLikeOtt,pOtt,iterOtt]=EM(X,K,Rndstart)
%
% I = number of objects
% J = number of variables
% K = number of clusters of the partition
%
% maxiter=max number of iterations
%
% Rndstart = number of initial random starts
%
maxiter=1000; % 
I = size(X,1);
J = size(X,2);
epsilon=0.0000001;

%seed=1;             %you can restart with a fixed seed of the random numbers
%rand('state',seed)    %
%
for loop=1:Rndstart
    
 % initial partition U0 is given
  U0=randPU(I,K);
  %
   su=sum(U0);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given U, compute initial mean vector (compute centroids) %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Xmean0 = diag(1./su)*U0'*X;
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given mean vector and U0, compute initial Var-Covar Matrix  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   Sig0=zeros(J);
   for k=1:K
      for i=1:I,
         Sig0=Sig0+U0(i,k)*((X(i,:)-Xmean0(k,:))'*(X(i,:)-Xmean0(k,:)));
      end
   end
   Sig0=(1/I)*Sig0;
   Sig0m1=inv(Sig0);
   detSig0=det(Sig0);
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given U0, compute the initial prior probabilities   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   p0=(1/I)*sum(U0)';
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given components computes posterior probabilities   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   U0=zeros(I,K);
   for i=1:I
            for k=1:K
                U0(i,k) = U0(i,k) + ((2*pi).^-(J./2))*(detSig0.^-0.5)*exp(-0.5*((X(i,:)-Xmean0(k,:))*Sig0m1*(X(i,:)-Xmean0(k,:))'))*p0(k);
            end
        end
   U0=diag(1./sum(U0'))*U0;
   su=sum(U0);
   Xmean0 = diag(1./su)*U0'*X;
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % compute intitial Log-Likelihhod    %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %
   %
   LLike0=0;
   for i=1:I
      for k=1:K
                LLike0=LLike0 + U0(i,k)*((-(J./2)*log(2*pi))-(0.5*log(detSig0))-0.5*((X(i,:)-Xmean0(k,:))*Sig0m1*(X(i,:)-Xmean0(k,:))')+ log(p0(k))) - U0(i,k)*log(U0(i,k));
      end
   end

   for iter=1:maxiter
       %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % EXPECTATION Compute Posterior Probabilities                       %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %
        Sig0m1=inv(Sig0);
        detSig0=det(Sig0);
        U=zeros(I,K);
        for i=1:I
            for k=1:K
                U(i,k) = U(i,k) + ((2*pi).^-(J./2))*(detSig0.^-0.5)*exp(-0.5*((X(i,:)-Xmean0(k,:))*Sig0m1*(X(i,:)-Xmean0(k,:))'))*p0(k);
            end
        end
        su=sum(U);
        while sum(su==0)>0,
            'entra'
            [m,p1]=min(su);
            [m,p2]=max(su);
            ind=find(U(:,p2));
            ind=ind(1:floor(su(p2)/2));
            U(ind,p1)=1;
            U(ind,p2)=0;
            su=sum(U);
        end
        U=diag(1./sum(U'))*U;
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAXIMIZATION                                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % given U, compute Xmean (compute centroids) ans Variance-Covariance matrix %
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
        Sig=(1/I)*Sig;
        Sigm1=inv(Sig);
        detSig=det(Sig);
        

%      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % given U, compute the prior probabilities   %   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   % 
        p=(1/I)*sum(U)';
   %
   % 
   % compute Log-Likelihood function
   %
        LLike=0;
        for i=1:I
            for k=1:K
                 LLike=LLike + U(i,k)*((-(J./2)*log(2*pi))-(0.5*log(detSig))-0.5*((X(i,:)-Xmean(k,:))*Sigm1*(X(i,:)-Xmean(k,:))')+ log(p(k))) - U(i,k)*log(U(i,k));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%   
        %       stopping rule     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%% 
        dif=LLike-LLike0;
        if dif < 0 & abs(dif)> epsilon 
            disp(sprintf('EM: Log-Likelihood decreases iter=%g, LLike0=%g, dif=%g', iter, LLike0, dif))
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
   disp(sprintf('EM: loop=%g, LLike=%g, iter=%g',loop,LLike,iter))
   if LLike > LLikeOtt
        UOtt=U;
        LLikeOtt=LLike;
        loopOtt=loop;
        iterOtt=iter;
        XmeanOtt=Xmean;
        pOtt=p;
   end
end
disp(sprintf('EM GaussianMM: loopOtt=%g, fOtt=%g, iter=%g',loopOtt,LLikeOtt,iter))