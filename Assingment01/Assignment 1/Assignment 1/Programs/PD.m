%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parsimonious Dendrogram %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Vichi Maurizio 
% version 20.10.2012
% 
% model:
%
% min||D-U*Db*U'-U*Dw*U'+diag(diag(U*Dw*U'))||^2
%
% subject to 
% 
% U binary and row stochastic
% Db ultrametric
%
function [UOtt,DbOtt, DwOtt, fOtt,iterOtt]=PD(D, K, Rndstart);

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



% total sum of squares
%
ts=sum(sum(D.^2));
%
% convergence tollerance
eps=0.000001;
% actual number of iterations 


Onesn=ones(n,1);
IK=eye(K);


for loop=1:Rndstart
    U=randPU(n,K);
    suu=sum(U)';
    ss=(suu==0);
    su=suu+ss;
    it=0;

    % ititial Dw 
    
    Dw=diag(diag(U'*D*U))*pinv((diag(1/suu)).^2-diag(1/suu));

    % initial Db
    Den=su*su';
    Db=(U'*D*U)./Den;
    Db=Db-diag(diag(Db));

    % Db ultrametric
    [Db]=upgma(Db,1);

    % inital value of the objective function 
    Q = U*Db*U'+U*Dw*U'-diag(diag(U*Dw*U'));
    fo=trace((D-Q)'*(D-Q)) / ts;
%
% Reiteration steps
%
    fdif=2*eps;
    while fdif > eps | it>=maxiter,
        it=it+1;
        %
        % update U
        %
        for i=1:n
            U(i,:)=IK(1,:);
            Q = U*Db*U'+U*Dw*U'-diag(diag(U*Dw*U'));
            miff=trace((D-Q)'*(D-Q)) / ts;
            posmin=1;
            for p=2:K
                U(i,:)=IK(p,:);
                Q = U*Db*U'+U*Dw*U'-diag(diag(U*Dw*U'));
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

   
        % update Db
        %
        su=sum(U)';
        ss=(su==0);
        su=su+ss;
    
        Den=su*su';
        Db=(U'*D*U)./Den;
        Db=Db-diag(diag(Db));
 
        [Db]=upgma(Db,1); 
        %
        %
        % update Dw
        %   
        %
        Dw=diag(diag(U'*D*U))*pinv((U'*U).^2-U'*U);
        %
        % check for convergence
        Q = U*Db*U'+U*Dw*U'-diag(diag(U*Dw*U'));
        f=trace((D-Q)'*(D-Q)) / ts;
        fdif=fo-f;
        fo=f;
    end
        disp(sprintf('Parsimonious Dendrogram: loop=%g, f=%g itr=%g',loop,f,it))
        if fdif<0
            sprintf('Parsimonious Dendrogram: error! fdif=%g negative: iteration=%g',fdif, it)
        end
         if loop==1
            UOtt=U;
            fOtt=f;
            loopOtt=1;
            iterOtt=1;
            DwOtt=Dw;
            DbOtt=Db;
        end
   if f < fOtt
        UOtt=U;
        fOtt=f;
        loopOtt=loop;
        iterOtt=it;
        DwOtt=Dw;
        DbOtt=Db;
   end
end
        disp(sprintf('Parsimonious Dendrogram (Final):  loopOtt=%g, fOtt=%g, iter=%g',loopOtt,fOtt,iterOtt))

        
