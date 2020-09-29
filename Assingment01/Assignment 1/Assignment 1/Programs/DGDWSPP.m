function [fwspp,awspp,bwspp,ag,bg,ari,Uwspp,Ug]=DGDWSPP(n,K, nsample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generator of Dissimilarities WSPP %
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
awspp  = zeros(nsample,1);
bwspp = zeros(nsample,1);
ag  = zeros(nsample,1);
bg = zeros(nsample,1);
ari = zeros(nsample,1);
Uwspp=cell(nsample,1);
Ug=cell(nsample,1);

for i=1:nsample
    % Generate dissimilarity data
    % generate partition
    U=randPU(n,K);
    a1=rand(1);
    a2=rand(1)*4+6;
    % generate normal error
    E = randn(n);
    E=E-diag(diag(E));
    E =(E+E')./2;
    D = a1*(ones(n)-U*U')+a2*((U*U')-eye(n))+E*7;
    %
    [UOwspp,aOwspp, bOwspp, fOwspp,iterOtt]=WSPP(D, K, 20);
    % solve lable switching problem
    C=zeros(K);
    for k=1:K
        for g=1:K
             C(k,g)=sum((U(:,k)-UOwspp(:,g)).^2 );
        end
    end
    [Pr,unfeas,DD] = bghungar(-C);
    %
    Ug{i}=U;
    Uwspp{i} = UOwspp(:,Pr);
    ag(i)=a1;
    bg(i)=a2;
    fwspp(i) = fOwspp;
    awspp(i) = aOwspp;
    bwspp(i) = bOwspp;
    ari(i)=mrand(Uwspp{i}'*U);

end

 