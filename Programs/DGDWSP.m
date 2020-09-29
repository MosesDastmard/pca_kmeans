function [fwspp,Dbwsp,Dwwsp,Dbg, Dwg, ari,Uwsp,Ug]=DGDWSP(n,K, nsample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generator of Dissimilarities WSP  %
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
fwspp = zeros(nsample,1);
ari = zeros(nsample,1);
Uwsp=cell(nsample,1);
Dbg=cell(nsample,1);
Dwg=cell(nsample,1);
Ug=cell(nsample,1);
Dbwsp=cell(nsample,1);
Dwwsp=cell(nsample,1);

for i=1:nsample
    % Generate dissimilarity data
    % generate partition
    U=randPU(n,K);
    Dw=diag(rand(K,1));
    Db=rand(K);
    Db=Db-diag(diag(Db));
    Db=(Db+Db')/.2;
    Db=Db*4+2*(ones(K)-eye(K));
    % generate normal error
    E = randn(n);
    E=E-diag(diag(E));
    E =(E+E')./2;
    D = U*Db*U'+U*Dw*U'-diag(diag(U*Dw*U'))+E*0.01;
    %
    [UOwsp,DbOwsp, DwOwsp, fOwsp]=WSP(D, K, 10);
    % solve lable switching problem
    C=zeros(K);
    for k=1:K
        for g=1:K
             C(k,g)=sum((U(:,k)-UOwsp(:,g)).^2 );
        end
    end
    [Pr,unfeas,DD] = bghungar(-C);    %
    Ug{i}=U;
    Dbg{i} = Db;
    Dwg{i} = Dw;
    Uwsp{i} = UOwsp(:,Pr);
    fwspp(i) = fOwsp;
    Dbwsp{i} = DbOwsp;
    Dwwsp{i} = DwOwsp;
    ari(i)=mrand(Uwsp{i}'*U);

end

 