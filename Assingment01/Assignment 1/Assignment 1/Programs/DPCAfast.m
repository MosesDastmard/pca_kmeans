%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disjoint Principal Component Analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maurizio Vichi
% October 2012

% problem maximize ||XBV||^2
% subject to
% A=BV where
% V'B'BV=I
% and V is binary and row stochastic
% B = diagonal matrix

%
function [Vdpca,Adpca, Ydpca fdpca, indpca]=DPCAfast(X,K,Rndstart)

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

% centring matrix
Jm=eye(n)-(1/n)*ones(n);

% compute var-covar matrix 
S=(1/n)*X'*Jm*X;

% Standardize data
Xs=Jm*X*diag(diag(S))^-0.5;


st=sum(sum(Xs.^2));

%um=ones(K,1);



for loop=1:Rndstart

    V=randPU(J,K);
    it=0;
    JJ=[1:J]';
    A=zeros(J,K);
    zJ=zeros(J,1);
    for g=1:K
            xx=V(:,g); 
            Jxx=[JJ(xx==1)];
            Sg=(1./n)*Xs(:,Jxx)'*Xs(:,Jxx);
            if sum(xx)>1
                [a,c]=eigs(Sg,1,'lm',opts);
                A(Jxx,g)=a;
            else
                A(Jxx,g)=1;
            end
    end
    A0=A;    
    Y=Xs*A;
    f0=trace((1./n)*Y'*Y);
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
                xx=V(:,g);           % classe attuale V
                xxx=V(:,posmax);     % classe vecchia V
                Jxx=[JJ(xx==1)];
                Jxxx=[JJ(xxx==1)];
                Sg=(1./n)*Xs(:,Jxx)'*Xs(:,Jxx);
                S0g=(1./n)*Xs(:,Jxxx)'*Xs(:,Jxxx);
                if sum(xx)>1
                     [a,c]=eigs(Sg,1,'lm',opts);
                      A(:,g)=zJ;
                      A(Jxx,g)=a;
                else
                      A(:,g)=zJ;
                      A(Jxx,g)=1;
                end
                if sum(xxx)>1
                      [aa,cc]=eigs(S0g,1,'lm',opts);
                      A(:,posmax)=zJ;  
                      A(Jxxx,posmax)=aa;
                else
                      A(:,posmax)=zJ;  
                      A(Jxxx,posmax)=1;
                end            
                Y=Xs*A;
                f=trace((1./n)*Y'*Y);
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
        f=trace((1./n)*Y'*Y);
        fdif = f-f0;
        
        if fdif > eps 
            f0=f; A0=A; 
        else
            break
        end
    end
  disp(sprintf('DPCA: Loop=%g, Explained variance=%g, iter=%g, fdif=%g',loop,f, it,fdif))   
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
[c,ic]=sort(varYdpca, 'descend');
Adpca=Adpca(:,ic);
Vdpca=Vdpca(:,ic);
Ydpca=Ydpca(:,ic); 

disp(sprintf('DPCA (Final): Explained variance=%g, loopdpca=%g, iter=%g, fdif=%g',fdpca, loopdpca, indpca,fdifo))

 