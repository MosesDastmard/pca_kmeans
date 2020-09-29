function [CIpF,PerOut] =StatBotStKm(XmOttd,UOttd,nus, pFd, PerIn) 
[G,Q,nrs]=size(XmOttd);
PerOut=zeros(G,Q,2);
%pValM=zeros(G,Q);

%MOttns=pValM;
CIpF=prctile(pFd,PerIn);
for k=1:G
    for j=1:Q
        PerOut(k,j,:)=prctile(XmOttd(k,j,:),PerIn);
%        MeOn(k,j)=mean(reshape(XmOttd(k,j,:),nrs,1));
%        VarOn(k,j)=cov(reshape(XmOttd(k,j,:),nrs,1),1);
        figure
        hist(reshape(XmOttd(k,j,:),nrs,1),30)
        figure
        normplot(reshape(XmOttd(k,j,:),nrs,1))
    end
end

%[MOttns, IX]=sort(MO);
%for j=1:Q
%    for k=1:G-1
%        tstI=MOttns(k,j)- (MOttns(k+1,j)-MOttns(k,j))./2;
%        tstS=MOttns(k,j)+ (MOttns(k+1,j)-MOttns(k,j))./2;
%        pValM(IX(k,j),j)=sum(MOtt(IX(k,j),j,:) >= tstS) ./ nrs;
%        pValM(IX(k,j),j)=pValM(IX(k,j),j)+sum(MOtt(IX(k,j),j,:) <= tstI) ./ nrs;
%        if k == G-1
%            tstI=tstS;
%            tstS=MOttns(k+1,j)+ (MOttns(k+1,j)-MOttns(k,j))./2;
%            pValM(IX(k+1,j),j)=sum(MOtt(IX(k+1,j),j,:) >= tstS) ./ nrs;
%            pValM(IX(k+1,j),j)=pValM(IX(k+1,j),j)+sum(MOtt(IX(k+1,j),j,:) <= tstI) ./ nrs;
%        end
%    end
%end
%for k=1:G
%    for j=1:Q
%        tst=MO(k,j);
%        %tst=(MO(k,j)-MeOn(k,j))./sqrt(VarOn(k,j));
%        if tst > 0
%            pValM(k,j)=sum(MOttn(k,j,:) >= tst) ./ nrs;
%        else
%            pValM(k,j)=sum(MOttn(k,j,:) <= tst) ./ nrs;
%       end
%    end
%end

%pVal = sum(dpFn*(1+f) >= pFO) ./ nrs; %
% fixed a significance level alfa=0.05 the critical avalue of the test is computed 
%critVal = prctile(dpFn,95);
%for the empirical distribution of H1, find the percentage of test
% statistic values that fall above the critical value of the test at 5%
% significance level
%power = sum(dpF*(1+f) >= critVal) ./ nrs;

