% generate data bivariate normal
% Vichi Metodologia Statistica Avanzata 2012-13

% Gaussian Bibariate 1 
%         Y  X
MU1 =    [0 0];
SIGMA1 = [1 0.5 ;
          0.5 1];
      
% Gaussian Bibariate 2
%         Y  X
MU2 =    [6 0];
SIGMA2 = [1 0 ;
          0 1];

% Gaussian Bibariate 3
%         Y  X
MU3 =    [3 5.2];
SIGMA3 = [1 -.4 ;
          -.4 1];
  
X = [mvnrnd(MU1,SIGMA1,500);    %(1/6)
     mvnrnd(MU2,SIGMA2,1500);   %(1/2)
     mvnrnd(MU3,SIGMA3,1000)];  %(1/3)

% plot data
scatter(X(:,1),X(:,2),10,'.')
hold on

pause

% compute AIC for differet number of components from 1 to 4 

AIC = zeros(1,4);
obj = cell(1,4);
for k = 1:4
    obj{k} = gmdistribution.fit(X,k);
    AIC(k)= obj{k}.AIC;
end

% the optimal number of components is given for the mimimum AIC  
[minAIC,numComponents] = min(AIC);
numComponents

% fit the gmm model for the optimal number of components
obj = gmdistribution.fit(X,numComponents);
h = ezcontour(@(x,y)pdf(obj,[x y]),[-4 10],[-4 10]);

% estimation of the mean vector var-covar matrices
ComponentMeans = obj.mu
ComponentCovariances = obj.Sigma
MixtureProportions = obj.PComponents


pause
hold off
figure

% surface of the GMM distribution
ezsurf(@(x,y)pdf(obj,[x y]),[-4 10],[-4 10])

pause

% plot of the clusters 
figure
scatter(X(:,1),X(:,2),10,'.')
hold on

idx = cluster(obj,X);
cluster1 = X(idx == 1,:);
cluster2 = X(idx == 2,:);
cluster3 = X(idx == 3,:);

delete(h)
h1 = scatter(cluster1(:,1),cluster1(:,2),10,'r.');
h2 = scatter(cluster2(:,1),cluster2(:,2),10,'g.');
h3 = scatter(cluster3(:,1),cluster3(:,2),10,'b.');

legend([h1 h2 h3],'Cluster 1','Cluster 2','Cluster 3','Location','NW')