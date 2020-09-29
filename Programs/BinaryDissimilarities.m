% Dissimilarities for binary variables 25.1.2008
% 
function [D] = BinaryDissimilarities(X, tp)
% x1 = binary codification of a binary variable 1
% x2 = binary codification of a binary variable 2
% ....
% xJ = binary codification of a binary variable J
% X=[x1, x2,...,xJ] binary data matrix (I x 2J)
%
% dissimilarity between objects xi and xl
[I,J]=size(X);
J=J/2;
D=zeros(J);
for i=1:I-1
    for l=i+1:I
                Cj  = reshape(X(i,:),2,J)';
                Ck  = reshape(X(l,:),2,J)';
                Cjk = Cj'*Ck;
                %pause
                if tp==1 
                    % Jaccard-Needham (c10 + c01)/(c11+c10+c01)
                    dil = (Cjk(2,1)+Cjk(1,2))/(Cjk(2,2)+Cjk(2,1)+Cjk(1,2));
                elseif tp==2
                    % Simple matching (c10 + c01)/ J
                    dil= (Cjk(2,1)+Cjk(1,2))/J;
                elseif tp==3
                    % Dice (c10 + c01)/(2*c11+c10+c01)
                    dil= (Cjk(2,1)+Cjk(1,2))/(2*Cjk(2,2)+Cjk(2,1)+Cjk(1,2));
                elseif tp==4
                    % Correlation coefficient [0.5 -(c11*c00-c10*c01)/((c10+c11)*(c01+c00)*(c11+c01)*(c00+c10))^0.5]
                    dil= 0.5 - (Cjk(2,2)*Cjk(1,1)- Cjk(2,1)*Cjk(1,2)) / ((Cjk(2,1)+Cjk(2,2))*(Cjk(1,2)+Cjk(1,1))*(Cjk(2,2)+Cjk(1,2))*(Cjk(1,1)+Cjk(2,1))).^0.5;
                elseif tp==5
                    % Yule (c10*c01)/(c11*c00+c10*c01)
                    dil = (Cjk(2,1)*Cjk(1,2))/(Cjk(2,2)*Cjk(1,1)+ Cjk(2,1)*Cjk(1,2));
                elseif tp==6
                    % Euclidean distance (C10+c01).^0.5
                    dil = (Cjk(2,1)+Cjk(1,2)).^0.5;
                elseif tp==7
                    % Variance (c10+c01)/(4*J)
                    dil = ((Cjk(1,2)+Cjk(2,1)))/(4*J);
                elseif tp==8
                    % Lance and Williams (c10+c01)/(2c00+c10+c01)
                    dil = ((Cjk(2,1)+Cjk(1,2)))/ (2*Cjk(1,1)+Cjk(1,2)+Cjk(2,1));
                elseif tp==9
                    % Rogers & Tanmoto (2*c10+2*c01) / (c00+c11+2*c10+2*c01)
                    dil = (2*Cjk(1,2)+2*Cjk(2,1))/ (Cjk(1,1)+Cjk(2,2)+2*Cjk(1,2)+2*Cjk(2,1));
                elseif tp==10
                    % Russell & Rao 1 - c11/J
                    dil = 1 - Cjk(2,2)/J;
                elseif tp == 11
                    % Kulzinsky (c10+c01-c11+J) /(c10+c01+J)
                    dil = (Cjk(1,2)+Cjk(2,1)-Cjk(2,2)+J) / (Cjk(1,2)+Cjk(2,1)+J);
                end
                D(i,l) = dil;
                D(l,i) = dil;
                
       end           
end
                
                    



%D=zeros(I);
%for i=1:I
%    for l=1:I

%        D(i,l) = BinDis(X,type);
%    end
%end


    

