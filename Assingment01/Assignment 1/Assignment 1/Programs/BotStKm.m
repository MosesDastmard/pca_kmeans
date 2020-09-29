%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap for Kmeans %
%%%%%%%%%%%%%%%%%%%%%%%%

% Vichi Maurizio
% 
% n = number of objects
% J = number of variables
% nrs = number of resampling

function [XmOttd, UOttd, nus, pFd,Uf]=BotStKm(X,U, nrs)
%  
%
[n,J]=size(X);
[n,G]=size(U);
%
% UOttd = Membership matrix bootstrap distribution
% XmOttd = Centroid matrix distribution
% nus = number of time each unit has been sampled in nrs resampling
% Uf = frequency of a unit classified in a cluster
% pFd = Pseudo F distribution
%
UOttd=zeros(n,G, nrs);
XmOttd=zeros(G,J,nrs);
nus=zeros(n,1);
Uf=zeros(n,G);
pFd=zeros(nrs,1);
su=sum(U);
Xm=diag(1./su)*U'*X;
%
% compute Pseudo F for the clustering of K-means
%
pF=(trace((U*Xm)'*(U*Xm))./(G-1))./(trace((X-U*Xm)'*(X-U*Xm))./(n-G))

% null hypothesis for testing absence of clustering structure

mu=mean(X);
SIGMA=cov(X,1);
DD=diag(diag(SIGMA).^0.5);
MR=corr(X);
m0=zeros(J,1);
for i=1:n
    Xin(i,:)=mvnrnd(m0,MR,1);
    while sum(find(Xin(i,:)>2 | Xin(i,:)<-2)) > 0
        Xin(i,:)=mvnrnd(m0,MR,1);
    end
end
Xin=Xin+ones(n,1)*mu;
Xin=Xin*DD;    
% Optimal K-means partitioning for null hypothesis
%
[loopOin,UOin,fOin,iterOin]=kmeansV(Xin,G,20);

    su=sum(UOin);
    Xmin=diag(1./su)*UOin'*Xin;
C=zeros(G);
    for k=1:G
        for g=1:G
%            C(k,g)=sum((MO(k,:)-MOn(g,:)).^2 );
             C(k,g)=sum((U(:,k)-UOin(:,g)).^2 );
        end
    end
    [Pr,unfeas,D] = BGHUNGAR(-C);
    
    Xmin(:,:)=Xmin(Pr,:);   

for i=1:nrs
%    i
% --------------- stratified random sampling without replacement -----
    f=0.8;
    dclas=U'*ones(n,1);
    dclasr=ceil(f*dclas);
    nr=dclasr'*ones(G,1);
    Samp=zeros(n,1);
    for g=1:G
        Str_g = find(U(:,g) == 1);
        PStr_g=Str_g(randperm(dclas(g)));
        Samp(PStr_g(1:dclasr(g)))= 1;
    end
    nus=nus+Samp;
    Xs=X(find(Samp==1),:);
% ---------------- random sample with replacement ---------------------
%   RSR = (ceil(n*rand(n,1)));
%   Xs=X(RSR,:);
% ---------------------------------------------------------------------
    [lpOs,UOs,fOs,iterOs]=kmeansV(Xs,G,20);

%    [result,UOt,VOt,MOt]=simuDKMA(Xs,G,Q,10);
%   [lOt,UOt,fOt,iterOtt]=kmeansVICHI(Xs,G,30);
    su=sum(UOs);
    Xms=diag(1./su)*UOs'*Xs;
    pFd(i)=(trace((UOs*Xms)'*(UOs*Xms))./(G-1))./(trace((Xs-UOs*Xms)'*(Xs-UOs*Xms))./(nr-G));
    U1=zeros(n,G);
    U1(find(Samp==1),:)=UOs;
    C=zeros(G);
        for k=1:G
            for g=1:G
%               C(k,g)=sum(( MO(k,:)-MOt(g,:) ).^2 );
                C(k,g)=sum(( U(:,k)-U1(:,g) ).^2 );
            end
        end
    % find the best matching for lable switching of units    
    [Pr,unfeas,D] = BGHUNGAR(-C);
    UOttd(find(Samp==1),:,i)=UOs(:,Pr);

    XmOttd(:,:,i)=Xms(Pr,:);      
    Uf=Uf+UOttd(:,:,i);
    
% ---------------------------------------------------------------------
%   null model for the sample
%
%    for ii=1:nr
%        Xsin(ii,:)=mvnrnd(m0,MR,1);
%        while sum(find(Xsin(ii,:)>2 | Xsin(ii,:)<-2)) > 0
%            Xsin(ii,:)=mvnrnd(m0,MR,1);
%        end
%    end
%    Xsin=Xsin+ones(nr,1)*mu;
%    Xsin=Xsin*DD;    

% ---------------------------------------------------------------------

%    [lpOsin,UOsin,fOsin,itOsin]=kmeansV(Xs,G,20)
%    su=sum(UOsin);
%    Xmsin=diag(1./su)*UOsin'*Xsin;
    
%    [result,UOtn,VOtn,MOtn]=simuDKMA(Xsn,G,Q,5);
%   [lOn,UOn,fOn,iterOttn]=kmeansVICHI(Xsn,G,10);
%    pFsd(i)=(trace((UOsin*Xmsim)'*(UOsin*Xmsim))./(G-1))./(trace((Xsin-UOsin*Xmsim)'*(Xsin-UOsin*Xmsim))./(nr-G));
%     U1=zeros(nr,G);
%     U1(find(Samp==1),:)=UOsin;
%    C=zeros(G);
%        for k=1:G
%            for g=1:G
%%%                C(k,g)=sum(( MOn(k,:)-MOtn(g,:) ).^2 );
%                C(k,g)=sum(( UOn(k,:)-U1(g,:) ).^2 );
%            end
%        end
%
%    % find the best matching for lable switching of units    
%    [Pr,unfeas,D] = BGHUNGAR(-C);
%    UOttn(find(Samp==1),:,i)=UOtn(:,Pr);
  %  MOttn(:,:,i)=MOtn(Pr,:);      
end
for i=1:n;  Uf(i,:)=Uf(i,:)./nus(i);end

%hist(dfO,30)
%figure
%hist(dpFn,30)
%figure
%hist(dpF,30)