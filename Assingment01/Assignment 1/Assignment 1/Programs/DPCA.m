%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disjoint Principal Component Analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maurizio Vichi
% October 2012 revised June 2013

% problem maximize ||XBV||^2
% subject to
% A=BV where
% V'B'BV=I
% and V is binary and row stochastic
% B = diagonal matrix

%
function [Vdpca,Adpca, Ydpca fdpca, indpca]=DPCA(X,K,Rndstart)

% X (n X J) data matrix
% V (J x K) membership matrix for clustering variables
% n = numenr of objects
% J = number of variables
% K = numberr of 



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

% compute var-covar matrix 
S=cov(Xs,1);

st=(1./n)*sum(sum(Xs.^2));


for loop=1:Rndstart

    V=randPU(J,K);
    it=0;
    JJ=[1:J]';
    A=zeros(J,K);
    zJ=zeros(J,1);
    for g=1:K
            ibCg=V(:,g); 
            JCg=[JJ(ibCg==1)];
            Sg=S(JCg,JCg);
            if sum(ibCg)>1
                [a,c]=eigs(Sg,1,'lm',opts);
                A(JCg,g)=a;
            else
                A(JCg,g)=1;
            end
    end
    A0=A;
    Y=Xs*A;
    f0=trace((1./n)*(Y'*Y));
    fmax=0;

    % iteration phase
    fdif=2*eps;
    while fdif > eps | it>=maxiter,
        it=it+1;
           
        % update V and A
        for j=1:J
            posmax=JJ(V(j,:)==1);
            for g=1:K
                V(j,:)=VC(g,:);
                ibCg=V(:,g);           % classe attuale V
                ibCpm=V(:,posmax);     % classe vecchia V
                JCg=[JJ(ibCg==1)];
                JCpm=[JJ(ibCpm==1)];
                Sg=S(JCg,JCg);
                S0g=S(JCpm,JCpm);
                if sum(ibCg)>1
                     [a,c]=eigs(Sg,1,'lm',opts);
                      A(:,g)=zJ;
                      A(JCg,g)=a;
                else
                      A(:,g)=zJ;
                      A(JCg,g)=1;
                end
                if sum(ibCpm)>1
                      [aa,cc]=eigs(S0g,1,'lm',opts);
                      if sum(aa)<0
                        aa=-aa;
                      end

                      A(:,posmax)=zJ;  
                      A(JCpm,posmax)=aa;
                else
                      A(:,posmax)=zJ;  
                      A(JCpm,posmax)=1;
                end            
                Y=Xs*A;
                f=trace((1./n)*(Y'*Y));
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
        
        Y=Xs*A; 
        f=trace((1./n)*(Y'*Y));
        fdif = f-f0;
        
        if fdif > eps 
            f0=f; fmax=f0;A0=A; 
        else
            break
        end
    end
  disp(sprintf('DPCA: Loop=%g, Explained variance=%g, iter=%g, fdif=%g',loop,f./st*100, it,fdif))   
       if loop==1
            Vdpca=V;
            Adpca=A;
            Ydpca=Xs*Adpca;
            fdpca=f;
            loopdpca=1;
            indpca=it;
            fdifo=fdif;
        end
   if f > fdpca
       Vdpca=V;
       fdpca=f;
       Adpca=A;
       Ydpca=Xs*Adpca;
       loopdpca=loop;
       indpca=it;
       fdifo=fdif;
   end
end
% sort the final solution in descend order of variance
varYdpca=var(Ydpca,1);
[~,ic]=sort(varYdpca, 'descend');
Adpca=Adpca(:,ic);
Vdpca=Vdpca(:,ic);
Ydpca=Ydpca(:,ic); 

disp(sprintf('DPCA (Final): Explained variance=%g, loopdpca=%g, iter=%g, fdif=%g',fdpca./st*100, loopdpca, indpca,fdifo))

 