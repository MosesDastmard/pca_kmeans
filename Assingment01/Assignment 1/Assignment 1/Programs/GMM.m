
%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Gaussian Mixture Model %
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% finds the best fit of the Gaussian Mixture Model, choosing the
% number of cluster by AIC 
function [IC,CM,CC,MP] = GMM(X);
% mnc = max number of components
% nc  = number of compontns
% ComponentMeans = centroid matrix (means vectors of the components)
% ComponentCovariances = matrices of covariances for each component
% MixtureProportions = vector of the mixture proportions
%
% compute AIC for differet number of components from 1 to 6 
mnc = 6;

AIC = zeros(1,mnc);
obj = cell(1,mnc);
for k = 1:mnc
    obj{k} = gmdistribution.fit(X,k, 'replicates',5);
    AIC(k)= obj{k}.AIC;
end

[minAIC,nc] = min(AIC);
nc
% fit the gmm model for the optimal number of components

% results of estimation of the mean vector var-covar matrices
CM = obj{nc}.mu;
CC = obj{nc}.Sigma;
MP = obj{nc}.PComponents;
% index of the clusters 
IC = cluster(obj{nc},X);
end

