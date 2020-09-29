%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disjoint Probabilistic Principal Component Analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maurizio Vichi
% May 2013
%
%
% X = YV'B + E 
% with covariance structure
% Sx = BVV'B + s2 I
%
% subject to
% A=BV where
% V'B'BV = I

% and V is binary and row stochastic
% s2 = variance of the error
% B = (J x J) diagonal matrix
% C = (K x K) diagonal matrix

% Maximum likelihood estimation by EM
%
%
% Maximize Likelihood -n/2 ( J*ln(2Pi) + ln(|Sx|) + tr(Sx^-1 S) )

%
function [Vdppca,Adppca, s2dppca, Ydppca fdppca, indppca]=DPPCA(X,K,Rndstart)

% X (n X J) data matrix
% S = (1/n)*X'*Jm*X; observed sample covariance matrix 
% Jm centering matrix
% V (J x K) membership matrix for clustering variables
% n = numenr of objects
% J = number of variables
% K = numberr of classes of variables 

%
% initialization
%
maxiter=100;
% convergence tollerance
eps=0.0000000001;

[n,J]=size(X);
VC=eye(K);
opts.disp=0;




% Standardize data
Xs=zscore(X,1);

% compute the observed var-covar matrix 
S=cov(Xs,1);

for loop=1:Rndstart

    V=randPU(J,K);
    rnV=sqrt(1./sum(V));
    it=0;
    % intial s2
    [~,L]=eig(S);
    dL=diag(L);
    s2=(1./(J-K))*sum(dL(K+1:J));   
    
    % initial A (A0)
    %
    JJ=[1:J]';
    A0=zeros(J,K);
    A=A0;
    zJ=zeros(J,1);
     
    for g=1:K
        ibCg=V(:,g); 
        JCg=[JJ(ibCg==1)];
        Sg=S(JCg,JCg);
        if sum(ibCg)>1
            [a,c]=eigs(Sg,1,'lm',opts);
            A0(JCg,g)=a*c.^0.5-s2^0.5;
        else
            A0(JCg,g)=1;
        end
    end
    % intial Sx
    Sx = A0*A0'+s2*eye(J);  
    f0 = (-n/2)*( J*log(pi) + log(det(Sx)) + trace(Sx\S) );
    fmax=f0;

    % iteration phase
    fdif=2*eps;
    while fdif > eps || it>=maxiter,
        it=it+1;
           
        % update V and A
        for j=1:J
            posmax=JJ(V(j,:)==1);
            for g=1:K
                V(j,:)=VC(g,:);
                ibCg=V(:,g);           % new class of V
                ibCpm=V(:,posmax);     % old class of V
                JCg=[JJ(ibCg==1)];
                JCpm=[JJ(ibCpm==1)];
                Sg=S(JCg,JCg);
                S0g=S(JCpm,JCpm);
                if sum(ibCg)>1
                     [a,c]=eigs(Sg,1,'lm',opts);
                      A(:,g)=zJ;                     
                      if sum(a)<0
                          a=-a;
                      end
                      A(JCg,g)=a*(c-s2).^0.5;
                else
                      A(:,g)=zJ;
                      A(JCg,g)=(1-s2).^0.5;
                end
                if sum(ibCpm)>1
                     [aa,cc]=eigs(S0g,1,'lm',opts);
                     if sum(aa)<0
                        aa=-aa;
                     end

                     A(:,posmax)=zJ;  
                     A(JCpm,posmax)=aa*(cc-s2).^0.5;
                else
                      A(:,posmax)=zJ;  
                      A(JCpm,posmax)=(1-s2).^0.5;
                end            
                Sx = A*A' + s2*eye(J);
                f = (-n/2)*( J*log(pi) + log(det(Sx)) + trace(Sx\S) );
                if f > fmax
                    fmax=f;
                    posmax=g;
                    A0=A;
                else
                    A=A0;
                end
            end
            V(j,:)=VC(posmax,:);
        end
        
        % update s2
        
        s2=0;
        for g=1:K
            ibCg=V(:,g); 
            JCg=[JJ(ibCg==1)];
            Sg=S(JCg,JCg);
            [U,L]=eig(Sg);
            [l,il]=sort(diag(L), 'descend');
            L=diag(l);
            U=U(:,il);
            s2 = s2+sum(l(2:size(Sg,1)));
        end
        s2=s2./(J-K);
        Sx = A*A' + s2*eye(J);
        f=(-n/2)*( J*log(pi) + log(det(Sx)) + trace(Sx\S) );

        fdif = f-f0;
        
        if fdif > eps 
            f0=f; fmax=f0; A0=A;
        else
            break
        end
    end
  disp(sprintf('DPPCA: Loop=%g, Likelihood function=%g, iter=%g, fdif=%g',loop,f, it,fdif))   
       if loop==1
            Vdppca=V;
            Adppca=A;
            s2dppca=s2;
            fdppca=f;
            loopdppca=1;
            indppca=it;
            fdifo=fdif;
        end
   if f > fdppca
       Vdppca=V;
       fdppca=f;
       Adppca=A;
       s2dppca=s2;
       loopdppca=loop;
       indppca=it;
       fdifo=fdif;
   end
end

% compute factors
%
% normalize factor loading matrix
Adppca=Adppca*inv(Adppca'*Adppca)^0.5;
Ydppca=Xs*Adppca;
% compute orthgonal factors
[C,L]=eig(cov(Ydppca,1));
Ydppca=Ydppca*C;

% sort the final solution in descend order of variance
varYdppca=var(Ydppca,1);
[~,ic]=sort(varYdppca, 'descend');
Adppca=Adppca(:,ic);
Vdppca=Vdppca(:,ic);
Ydppca=Ydppca(:,ic); 

disp(sprintf('DPPCA (Final): Likelihood function=%g, loopdpca=%g, iter=%g, fdif=%g',fdppca, loopdppca, indppca,fdifo))

 