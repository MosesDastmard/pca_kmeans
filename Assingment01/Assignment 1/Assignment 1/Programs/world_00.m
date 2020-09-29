function [B,T,f]=varim(A);
% produces varimax rotated version of A and rotation matrix T
% Based on nevels 1986; see also pascal source
%
% input: A    (matrix to be rotated)
% output: B (rotated A)
%	  T (rotation matrix)
%	  f (varimax function)


conv=.000001;
[m,r]=size(A);
%T=eye(r);
T=zeros(r,r);
for i=1:r,T(i,i)=1;end;
B=A;

f=ssq((A.*A)-ones(m,1)*sum(A.*A)/m);
fold=f-2*conv*f;
if f==0, fold=-conv; end;
iter=0;

while f-fold>f*conv 
  fold=f;iter=iter+1;
  for i=1:r
    for j=i+1:r
      x=B(:,i);y=B(:,j);
      xx=T(:,i);yy=T(:,j);
      u=x.^2-y.^2;v=2*x.*y;
      u=u-ones(m,1)*sum(u)/m;v=v-ones(m,1)*sum(v)/m;  
      a=2*sum(u.*v);b=sum(u.^2)-sum(v.^2);c=(a^2+b^2)^.5;
      if a>=0; sign=1; end;
      if a<0; sign=-1; end;  
      if c<.00000000001 
%        disp(' No rotation anymore');
        cos=1;sin=0;
      end;
      if c>=.00000000001 
        vvv=-sign*((b+c)/(2*c))^.5;
        sin=(.5-.5*vvv)^.5;cos=(.5+.5*vvv)^.5;
      end;
      v=cos*x-sin*y;w=cos*y+sin*x;
      vv=cos*xx-sin*yy;ww=cos*yy+sin*xx;
      if vvv>=0       % prevent permutation of columns
        B(:,i)=v;B(:,j)=w;T(:,i)=vv;T(:,j)=ww;
      end;
      if vvv<0
        B(:,j)=v;B(:,i)=w;T(:,j)=vv;T(:,i)=ww;
      end;
    end;
  end;  
  f=ssq((B.*B)-ones(m,1)*sum(B.*B)/m)  ;
%mbscalar(cos);
%mbscalar(sin);
%mbscalar(a);
%mbscalar(b);
%mbscalar(c);
%mbscalar(f);
%mbvector(w);
%mbvector(ww);
%mbvector(x);
%mbvector(xx);
%mbvector(y);
%mbvector(yy);
%mbvector(u);
%mbvector(v);
%mbvector(vv);


end;

function t = ssq(a)
%SSQ	SSQ(A) is the sum of squares of the elements of matrix A.
% Written by Henk A.L. Kiers, University of Groningen (Last update: September 1, 1999)

t = sum(sum(a.^2));