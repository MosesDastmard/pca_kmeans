function [ari,Ug,Umlc,Umlmc,Uem]=XGMBC(n,Pr, nsample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X Generator for model based clustering  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% n = number of objects
% K = number of clusters
% nsample = number of samples
%
% Rndstart number of random starts
%
%
K=size(Pr,1);
fwspp = zeros(nsample,1);
ari = zeros(nsample,3);
Ug=cell(nsample,1);
Umlc=cell(nsample,1);
Umlmc=cell(nsample,1);
Uem=cell(nsample,1);


for i=1:nsample
    'campione', i
    % Generate dissimilarity data
    % generate partition
    [X,U,u]=MixSampling(n, Pr);
    K=size(Pr,1);
    
    % Compute model MLC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [loopOtt,UOmlc,LLikeOtt,iterOtt]=MLC(X,K,10);
    %
    % solve lable switching problem
    C=zeros(K);
    for k=1:K
        for g=1:K
             C(k,g)=sum((U(:,k)-UOmlc(:,g)).^2 );
        end
    end
    [Pm,unfeas,DD] = bghungar(-C);    %
    Ug{i}=U;
    Umlc{i} = UOmlc(:,Pm);
    ari(i,1)=mrand(Umlc{i}'*U);
    % Finish model MLC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute model MLMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        [loopOtt,UOmlmc,LLikeOtt,pOtt,iterOtt]=MLMC(X,K,10);
    %
    % solve lable switching problem
    C=zeros(K);
    for k=1:K
        for g=1:K
             C(k,g)=sum((U(:,k)-UOmlmc(:,g)).^2 );
        end
    end
    [Pm,unfeas,DD] = bghungar(-C);    %
    Umlmc{i} = UOmlmc(:,Pm);
    ari(i,2)=mrand(Umlmc{i}'*U);
    % Finish model MLMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [loopOtt,UOemf,LLikeOtt,pOtt,iterOtt]=EM(X,K,5);
    %
    % compute hard UOem
    [mU,inmuem]=max(UOemf');
    Ik=eye(K);
    UOem=Ik(inmuem,:);
    % solve lable switching problem
    C=zeros(K);
    for k=1:K
        for g=1:K
             C(k,g)=sum((U(:,k)-UOem(:,g)).^2 );
        end
    end
    [Pm,unfeas,DD] = bghungar(-C);    %
    Uem{i} = UOem(:,Pm);
    ari(i,3)=mrand(Uem{i}'*U);
    % Finish model MLMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

 