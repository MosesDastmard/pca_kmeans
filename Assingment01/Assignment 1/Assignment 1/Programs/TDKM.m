% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       %   
%       Trimmed Double K-Means          %
%                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Maurizio Vichi 18.8.2013 
%
% X (n X J) data matrix
% V (J x Q) membership matrix for clustering variables
% U (n x k) membership matrix for clustering objects
%
% model X = DUYmV'B + E
% 
% problem: min||X-DUYmV'B||^2
% 
% being ||X||^2 = ||X-UYmV'B||^2 + ||UYmV'B||^2
%
% equivalent problem
%
% problem maximize ||UYmV'B||^2
% subject to
% V binary and row stochastic
% U binary and row stochastic
% B Diagonal matrix binary 
% D Diagonal matrix binary
function [Vtdkm,Utdkm,dtdkm, btdkm,ftdkm,intdkm]=TDKM(X, K, Q, alfaO,alfaV, Rndstart)


% 
% initialization
%
maxiter=100;
% convergence tollerance
eps=0.0000000001;

opts.disp=0;
VC=eye(Q);
[n,J]=size(X);

d=ones(n,1);
done=d;
dmin=zeros(n,1);
b=ones(J,1);
bone=b;
bmin=zeros(J,1);
D=eye(n);
B=eye(J);
% Standardize data
Xs=zscore(X,1);

% compute var-covar matrix 

S=cov(Xs,1);
st = sum(sum(Xs.^2));


for loop=1:Rndstart
    V=randPU(J,Q);
    U=randPU(n,K);
    % update centroid matrix Xmean 
    su=sum(U);
    sv=sum(V);
    Ym = diag(1./su)*U'*Xs*V*diag(1./sv);
    % compute na and Ja
    na=round(n*alfaO);
    Ja=round(J*alfaV);
    %
    % take the first na object as those with maximal distance from centroid
    nadmax=zeros(na,1);
    pnadmax=zeros(na,1);
    Ymv=Ym*V';
    for i=1:na
        nadmax(i)=sum((Xs(i,:)-Ymv(1,:)).^2);
        pnadmax(i)=i;
        for k=2:K
            dif=sum((Xs(i,:)-Ymv(k,:)).^2);
            if dif < nadmax(i)
                nadmax(i)=dif;
            end 
        end
    end
    % find the na objects with maximal distance
    [nadmin, pnadmin]=min(nadmax);
    for i=na+1:n
        mindif=sum((Xs(i,:)-Ymv(1,:)).^2);
        for k=2:K
            dif=sum((Xs(i,:)-Ymv(k,:)).^2);
            if dif < mindif
                mindif=dif;
            end 
        end
        if mindif > nadmin
            nadmax(pnadmin)=mindif;
            pnadmax(pnadmin)=i;
            [nadmin, pnadmin]=min(nadmax);
        end
    end
    %
    % update D
    %
    d=done;
    d(pnadmax)=0;
    D=diag(d);
    ind=find(d==1);
    %
    % take the first Ja variables as those with maximal distance from centroid
    Ymu=U*Ym;
    Jadmax=zeros(Ja,1);
    pJadmax=zeros(Ja,1);
    for j=1:Ja
        Jadmax(j)=sum((Xs(:,j)-Ymu(:,1)).^2);
        pJadmax(j)=j;
        for q=2:Q
            dif=sum((Xs(:,q)-Ymu(:,q)).^2);
            if dif < Jadmax(j)
                Jadmax(j)=dif;
            end 
        end
    end
    % find the Ja variables with maximal distance
    [Jadmin, pJadmin]=min(Jadmax);
    for j=Ja+1:J
        mindif=sum((Xs(:,j)-Ymu(:,1)).^2);
        posmin=j;
        for q=2:Q
            dif=sum((Xs(:,q)-Ymu(:,q)).^2);
            if dif < mindif
                mindif=dif;
            end 
        end
        if mindif > Jadmin
            Jadmax(pJadmin)=mindif;
            pJadmax(pJadmin)=posmin;
            [Jadmin, pJadmin]=min(Jadmax);
        end
    end
    %
    % update B
    %
    b=bone;
    b(pJadmax)=0;
    B=diag(b);
    inb=find(b==1);
    %
    % update centroid matrix Xmean 
    su=sum(D*U);
    sv=sum(B*V);
    Ym = diag(1./su)*U'*D*Xs*B*V*diag(1./sv);
    %
    % compute objective funcion
    BB=D*U*Ym*V'*B;
    %f0=trace(BB'*BB)/st;
    f0=trace(BB'*BB)/trace((D*Xs*B)'*(D*Xs*B));
    fmax=0;
    it=0;

%   iteration phase
    fdif=2*eps;
    while fdif > eps || it>=maxiter,
        it=it+1;

        % find objects and variables to be trimmed
        
        % find the na objects with maximal distance

        [nadmin, pnadmin]=min(nadmax);
       
        for i=1:n-na
            id=ind(i);
            midifu=sum((Xs(id,inb)-Ymv(1,inb)).^2);
            for k=2:K
                dif=sum((Xs(id,inb)-Ymv(k,inb)).^2);
                if dif < midifu
                    midifu=dif;
                end 
            end
            if midifu > nadmin
                nadmax(pnadmin)=midifu;
                pnadmax(pnadmin)=id;
                [nadmin, pnadmin]=min(nadmax);
            end
        end
        %
        % update D
        %
 

        d=done;
        d(pnadmax)=0;
        D=diag(d);
        ind=find(d==1);
        %
        % find the Ja variables with maximal distance
        %
        [Jadmin, pJadmin]=min(Jadmax);
        for j=1:J-Ja
            jb=inb(j);
            midifv=sum((Xs(ind,jb)-Ymu(ind,1)).^2);
            for q=2:Q
                dif=sum((Xs(ind,jb)-Ymu(ind,q)).^2);
                if dif < midifv
                    midifv=dif;
                end 
            end
            if midifv > Jadmin
                Jadmax(pJadmin)=midifv;
                pJadmax(pJadmin)=jb;
                [Jadmin, pJadmin]=min(Jadmax);
            end
        end
        %
        % update B
        %
        b=bone;
        b(pJadmax)=0;
        B=diag(b);
        inb=find(b==1);
        %
        %
        % given Ymean and V update U
        U=zeros(n,K);
        Ymv=Ym*V';
        for i=1:n-na
            id=ind(i);
            mdif=sum((Xs(id,inb)-Ymv(1,inb)).^2);
            posmin=1;
            for k=2:K
                dif=sum((Xs(id,inb)-Ymv(k,inb)).^2);
                if dif < mdif
                    mdif=dif;
                    posmin=k;
                end 
            end
            U(id,posmin)=1;
        end
     %
        su=sum(D*U);
        flgU=0;
        while sum(su==0)>0,
            %'entra U'
            flgU=1;
            [~,p1]=min(su);
            [~,p2]=max(su);
            indd=find(D*U(:,p2));
            indd=indd(1:floor(su(p2)/2));
            U(indd,p1)=1;
            U(indd,p2)=0;
            su=sum(D*U);
        end 

        % given U and V compute Ym (compute centroids)
        %
        su=sum(D*U);
        Ym = diag(1./su)*U'*D*Xs*B*V*diag(1./sv);
        %

        % given U and Ym update V
            V=zeros(J,Q);
            Ymu=U*Ym;
            for j=1:J-Ja
                jb=inb(j);
                mindif=sum((Xs(ind,jb)-Ymu(ind,1)).^2);
                posmin=1;
                for i=2:Q
                    dif=sum((Xs(ind,jb)-Ymu(ind,i)).^2);
                    if dif < mindif
                        mindif=dif;
                        posmin=i;
                    end 
                end
                V(jb,posmin)=1;
            end

            sv=sum(B*V);
            flgV=0;
            while sum(sv==0)>0,
                %'entra V'
                flgV=1;
                [~,p1]=min(sv);
                [~,p2]=max(sv);
                indd=find(B*V(:,p2));
                indd=indd(1:floor(sv(p2)/2));
                V(indd,p1)=1;
                V(indd,p2)=0;
                sv=sum(B*V);
            end 
            
       % given U and V updata Xm
       
        sv=sum(B*V);
        Ym = diag(1./su)*U'*Xs*V*diag(1./sv);
  
        BB=D*U*Ym*V'*B;
        %f=trace(BB'*BB)/st;
        f=trace(BB'*BB)/trace((D*Xs*B)'*(D*Xs*B));
        fdif = f-f0;
        if fdif > eps 
            f0=f;fmax=f0;
        else
            break
        end
    end
  disp(sprintf('TDKM: Loop=%g, Explained variance=%g, iter=%g, fdif=%g',loop,f.*100, it,fdif))   
       if loop==1
            Vtdkm=V;
            Utdkm=U;
            ftdkm=f;
            looptdkm=1;
            intdkm=it;
            dtdkm=d;
            btdkm=b;
            fdifo=fdif;
        end
   if f > ftdkm
       Vtdkm=V;
       Utdkm=U;
       Ymtdkm=Ym;
       ftdkm=f;
       dtdkm=d;
       btdkm=b;
       looptdkm=loop;
       intdkm=it;
       fdifo=fdif;
   end
end
% sort clusters of variables in descend order of variance

[~,ic]=sort(diag(Vtdkm'*Vtdkm), 'descend');
Vtdkm=Vtdkm(:,ic);

% sort clusters of objects in descending order of variance
%dwc=zeros(K,1);
[~,ic]=sort(diag(Utdkm'*Utdkm), 'descend');
Utdkm=Utdkm(:,ic);
disp(sprintf('TDKM (Final): Percentage Explained Variance=%g, looptdkm=%g, iter=%g, fdif=%g',ftdkm.*100, looptdkm, intdkm,fdifo))

 