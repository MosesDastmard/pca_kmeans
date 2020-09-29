%% 1) Using “MixSampling.m” function, generate 300 observations by a bivariate 
%     Gaussian mixture with a structure of three groups.

[X,U,u]=MixSampling(300, [0.33;0.33;0.34], 9, 0, 1.5);

%% 2) Apply the Well Structured Partition (WSP) and Well Structured Perfect 
%     Partition (WSPP) models on data generated in the step 1) fixing the 
%     number of clusters equal to 3. Compare the membership matrices obtained 
%     by the two models with the “real” memberships values. Comment the results.

D = squareform(pdist(X)); % distances matrix 


[UOtt_wspp,bOtt_wspp, aOtt_wspp, fOtt_wspp,iterOtt_wspp]=WSP(D, 3, 10); 
Cm = UOtt_wspp'*U; % confusion matrix 

%% 3) Generate 4 noise variables of 300 observations by a Gaussian 
%     distribution with m=0 and sigma=2. Add these 4 variables to the data 
%     generated in the step 1). Now, you have a new 300×6 data matrix.

Y = randn(300, 4)*2; %generate 4 noise variables
Z = [X Y]; % adding noise  

%% 4) Apply the k-means model on this new data set fixing the number of 
%     clusters equal to 3. Compare the membership matrix obtained by k-means 
%     with the “real” memberships values. Comment the results.

[loopOtt_km,UOtt_km,fOtt_km,iterOtt_km]=kmeansVICHI(Z,3,10);
Cm2 = UOtt_km'*U; % confusion matrix 

%% 5) “Using a dimensionality reduction technique, such as principal 
%     components analysis (PCA), to create new variables and then using 
%     cluster analysis, such as k-means, to form groups using these new 
%     variables”. This technique is called as “Tandem Analysis”. Use this 
%     approach on the new data set generated in the step 3) fixing the 
%     number of clusters equal to 3. Compare the membership matrix obtained 
%     by k-means with the “real” memberships values. Comment the results.


%% 6) Apply Reduced k-means (REDKM) on data generated in the step 3) fixing 
%     the number of clusters equal to 3. Compare the membership matrix 
%     obtained by the REDKM with the “real” memberships values. Comment the 
%     results underling the simultaneous approach advantages with respect to 
%     the sequential approach.

[Urkm,Arkm, Yrkm,frkm,inrkm]=REDKM(Z, 3, 2, 20);
Cm3 = Urkm'*U; % confusion matrix 

%% 7) Load MatLab workspace named “data02.mat”. Repeat the steps 5) and 6) 
%     on the ECSI data set. Define the partition obtained by both approaches 
%     and comment the results.

clear all; % clear all elements previously generated
load data02.mat; %loading workspace



