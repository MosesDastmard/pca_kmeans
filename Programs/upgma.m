function [U,A,B,levfus]=upgma(D,model)
% hierarchical classification 
% Vichi Maurizio
% UPGMA, SLINK, CLINK, WARD's
%
% nobj   = number of objects
% levfus = livels of fusion
% D      = matrix of dissimilarities(nobjxnobj)(input) 
% model  = 1 case: UPGMA
% model  = 2 case: SLINK
% model  = 3 case: CLINK
% model  = 4 case: WARD's
% A      = leaders A(istep) amalgamated at istep
% B      = leaders B(istep)
%
%
nobj  = size(D,1);
A  = zeros(nobj-1,1);
B  = zeros(nobj-1,1);
P  = zeros(nobj,1);
PP = zeros(nobj,1);
QQ = zeros(nobj,1);
for i=1:nobj,
%
% 'leader' of class to which object I belongs
%
	PP(i)=i;
%
% number of objects in class whose 'leader' is I
%
	QQ(i)=1;
end
% initial values 
%
istep = 1;
levfus=zeros(nobj-1,1);
%
% start the algorithm
%
nstep=nobj; % number of steps + 1 to which upgma has to be stopped 

for istep = 1:nstep-1,
   
  
% compute the minimum distance
%
   dmin=inf;
   for i=1:nobj-1,
      if PP(i) == i,
         for j=i+1:nobj,
            if PP(j) == j
               if D(i,j) < dmin 
                  ic=i;
                  jc=j;
                  dmin = D(i,j);
               end
            end
         end
      end
   end
       
  for j=1:nobj
    if PP(j) == jc
      PP(j)=ic;
    end
  end
%       write
%        'amalgamation', [istep, ic,jc]
%	PP
%
% save amalgamation step
%
% Amalgamation istep is between classes with
% leaders A(istep) and B(istep), and takes 
% place at level levfus(istep)
%
A(istep) = ic;
B(istep) = jc;	
levfus(istep)=dmin;

%       update distances 
%
  for i=1:nobj,
    if i ~= ic | PP(i) == i

      if model == 1, %case: +++ UPGMA +++
        ds=(QQ(ic).*D(ic,i)+QQ(jc).*D(jc,i))./(QQ(ic)+QQ(jc));
      end
      if model == 2, %case: +++ SLINK +++
        ds=min(D(ic,i),D(jc,i));
      end
      if model == 3, %case: +++ CLINK +++
        ds=max(D(ic,i),D(jc,i)); 
      end
      if model == 4, %case: +++ WARD's +++ 
        ds=(1./(QQ(ic)+QQ(jc)+QQ(i))).*((QQ(ic)+QQ(i)).*D(ic,i)+(QQ(i)+QQ(jc)).*D(jc,i)-QQ(i).*dmin);
      end
    D(ic,i) = ds;
    D(i,ic) = ds;
    end
  end
  QQ(ic)=QQ(ic)+QQ(jc);
  istep=istep+1;
end
% write results
% 'history of hierarchical classification'
% 'fusion level  classes aggregated'
% [levfus(:), A(:), B(:)]


% compute the ultrametric matrix U
% associated to matrix D according model i=1 or 2 or 3 or 4

for i = 1:nobj
  P(i) = i;
end
for k=1:nstep-1
  for i = A(k):nobj
    if P(i) == A(k), 
      for j = B(k):nobj,
        if P(j) == B(k),
          U(i,j) = levfus(k);
	       U(j,i) = levfus(k);
        end
      end
    end
  end
  for i=A(k):nobj
    if P(i) == B(k),
      P(i)=A(k);
    end
  end
end 
