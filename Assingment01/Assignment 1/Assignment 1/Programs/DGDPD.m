function [fPD,DbPD,DwPD,DbgPD,DwgPD,ari,UPD,Ug]=DGDPD(n,K, nsample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generator of Dissimilarities PD   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% n = number of objects
% K = number of clusters
% nsample = number of samples
%
% Rndstart number of random starts
%
%
fPD = zeros(nsample,1);
ari = zeros(nsample,1);
UPD=cell(nsample,1);
Ug=cell(nsample,1);
DbgPD=cell(nsample,1);
DwgPD=cell(nsample,1);
DbPD=cell(nsample,1);
DwPD=cell(nsample,1);

for i=1:nsample
    % Generate dissimilarity data
    % generate partition
    U=randPU(n,K);
    Dw=diag(rand(K,1));
    Db=rand(K);
    Db=Db-diag(diag(Db));
    Db=(Db+Db')/.2;
    Db=Db*4+2*(ones(K)-eye(K));
    
    % Db ultrametric
    [Db]=upgma(Db,1);
    
    % generate normal error
    E = randn(n);
    E=E-diag(diag(E));
    E =(E+E')./2;
    D = U*Db*U'+U*Dw*U'-diag(diag(U*Dw*U'))+E*0.01;
    %
    [UOwsp,DbOwsp, DwOwsp, fOwsp]=PD(D, K, 10);
    % solve lable switching problem
    C=zeros(K);
    for k=1:K
        for g=1:K
             C(k,g)=sum((U(:,k)-UOwsp(:,g)).^2 );
        end
    end
    [Pr,unfeas,DD] = bghungar(-C);    %
    Ug{i}=U;
    DbgPD{i} = Db;
    DwgPD{i} = Dw;
    UPD{i} = UOwsp(:,Pr);
    fPD(i) = fOwsp;
    DbPD{i} = DbOwsp;
    DwPD{i} = DwOwsp;
    ari(i)=mrand(UPD{i}'*U);

end

 