function [X,U,A,Ym,ElementiNonNulli]=randSparse2(N,J,K,Q,Eterogeneita,Isolamento)
% Si definisce un modello X=UYmA'+E
%epsilon definisce quanto siano "mischiati" i dati
E=randn(N,J)*Eterogeneita;
U=randPU(N,K);
Ym=rand(K,Q)*Isolamento;
A=calcoloA(Q,J);
X=U*Ym*A'+E;
ElementiNonNulli=zeros(1,Q);
for q=1:Q,
    ElementiNonNulli(q)=nnz(A(:,q));
end
end




%Funzione per avere una partizione casuale dei dati
function [U]=randPU(n,c)

% generates a random partition of n objects in c classes
%
% n = number of objects
% c = number of classes
%
U=zeros(n,c);
U(1:c,:)=eye(c);

U(c+1:n,1)=1;
for i=c+1:n
    U(i,[1:c])=U(i,randperm(c));
end
U(:,:)=U(randperm(n),:);
end

function A=calcoloA(Q,J)
%Calcolo casuale di A 
%a=zeros(K,Q);

A=randPU(J,Q);
v=A*[1:Q]';
[B,C]=sort(v);
A=A(C,:);
for j=1:J
    uno=rand(1)*(3);
    for q=1:Q
       if(rand(1)<0.1 && uno>0 && A(j,q)==0)
           A(j,q)=1;
           uno=uno-1;
       end           
     end
end

A=A.*rand(J,Q);
a=sum(A);
for j=1:Q,
    A(:,j)=A(:,j)/abs(a(j));          
end
A=A.^0.5;

end

