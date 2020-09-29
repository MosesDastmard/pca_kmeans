%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Well Structured Perfect Partition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Vichi Maurizio 
% version 20.10.2012
% 
% model:
%
% ||D-b(11'-UU')-a(UU'-I)||^2
%
% subject to 
% 
% U binary and row stochastic
% a  eterogeneity within clusters
% b isolation between clusters
%
function [UOtt,bOtt, aOtt, fOtt,iterOtt]=WSPP(D, K, Rndstart);

% n = number of objects
% K = number of clusters of the partition
%
%
% maxiter=max number of iterations
% D = matrix of (dis)similarities
% Rndstart = number of random starts
% 
% initialization
%
maxiter=100;
%
n=size(D,1);
eps=0.000001;


% total sum of squares
%
ts=sum(sum(D.^2));
%
% convergence tollerance
eps=0.0000000001;
% actual number of iterations 


Onesn=ones(n,1);
One=ones(n);
IK=eye(K);
In=eye(n);


for loop=1:Rndstart
    U=randPU(n,K);
    suu=sum(U)';
    ss=(suu==0);
    su=suu+ss;
    it=0;
    sumni2=sum(sum(U).^2);

    % ititial a 
    
     a=trace(U'*D*U)/(sumni2-n);
   

    % initial b

     b=(trace(Onesn'*D*Onesn)-trace(U'*D*U))/(n.^2 - sumni2);

    % inital value of the objective function 
    Q = b*(One-U*U')+a*(U*U'-In);
    fo=trace((D-Q)'*(D-Q)) / ts;
%
% Reiteration steps
%
    fdif=2*eps;
    while fdif > eps,
        it=it+1;
   % update U
   %
        for i=1:n
            U(i,:)=IK(1,:);
            Q = b*(One-U*U')+a*(U*U'-In);
            miff=trace((D-Q)'*(D-Q)) / ts;
            posmin=1;
            for p=2:K
                U(i,:)=IK(p,:);
                Q = b*(One-U*U')+a*(U*U'-In);
                ff=trace((D-Q)'*(D-Q)) / ts;
                if ff < miff
                    miff=ff;
                    posmin=p;
                end 
            end
            U(i,:)=IK(posmin,:);
        end
   %
        su=sum(U);
        while sum(su==0)>0,
            dw=U'*((D-U*pinv(U)*D*pinv(U)'*U').^2)*Onesn;
            [dw,su'];
            [m,p1]=min(su);
            [m,p2]=max(dw);
            ind=find(U(:,p2));
            ind1=ind(1:floor(su(p2)/2));
            U(ind1,p1)=1;
            U(ind1,p2)=0;
       %U1=kmeans2(X(ind,:),U(ind,[p1,p2]));
       %U(ind,[p1,p2])=U1;
            su=sum(U);
        end 

        sumni2=sum(sum(U).^2);

        % Updata  a 
    
        a=trace(U'*D*U)/(sumni2-n);
   

        % Updata b

         b=(trace(Onesn'*D*Onesn)-trace(U'*D*U))/(n.^2 - sumni2);
   
        %
        % check for convergence
        Q = b*(One-U*U')+a*(U*U'-In);
        f=trace((D-Q)'*(D-Q)) / ts;
        fdif=fo-f;
        fo=f;
    end
        disp(sprintf('Well Structured Perfect Partition: loop=%g, f=%g itr=%g',loop,f,it))
        if fdif<0
            sprintf('Well Structured Perfect Partition: error! fdif=%g negative: iteration=%g',fdif, it)
        end
         if loop==1
        UOtt=U;
        fOtt=f;
        loopOtt=1;
        iterOtt=1;
        aOtt=a;
        bOtt=b;
   end
   if f < fOtt
        UOtt=U;
        fOtt=f;
        loopOtt=loop;
        iterOtt=it;
        aOtt=a;
        bOtt=b;
   end
end
disp(sprintf('Well Structured Perfect Partition (Final):  loopOtt=%g, fOtt=%g, iter=%g',loopOtt,fOtt,iterOtt))

        
