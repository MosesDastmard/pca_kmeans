%% Step 1)
[X,U,u]=MixSampling(300,[.2;.3;.5], 9,0,1.5);
D = pdist(X);
D = squareform(D);
K=3;
Rndstart=10;
Q=2;
%% Step 2)
% [UOtt_P,DbOtt_P, DwOtt_P, fOtt_P,iterOtt_P]=WSP(D, K, Rndstart);
% Cm_P = UOtt_P'*U;
% [UOtt_PP,DbOtt_PP, DwOtt_PP, fOtt_PP,iterOtt_PP]=WSPP(D, K, Rndstart);
% Cm_PP = UOtt_PP'*U;
%% Step 3)
Y = randn(300,4)*2;
X_new = [X Y];


%% Step 4)
[loopOtt_KM,UOtt_KM,fOtt_KM,iterOtt_KM] = kmeansVICHI(X_new,K,Rndstart);
Cm_KM = UOtt_KM'*U;

% %% Step 5)
% [Coeff_PCA,Score_PCA,latent_PCA] = pca(X_new);
% Score_PCA = Score_PCA(:,1:6);
% Score_Rotated_PCA = zscore(Coeff_PCA);
% [loopOtt_PCAKM,UOtt_PCAKM,fOtt_PCAKM,iterOtt_PCAKM] = kmeansVICHI(Score_Rotated_PCA,K,Rndstart);
% Cm_PCAKM = UOtt_PCAKM'*U;

%% Step 6)
[Urkm_RKM,Arkm_RKM, Yrkm_RKM,frkm_RKM,inrkm_RKM] = REDKM(X_new, K, Q, Rndstart);
Cm_RKM = Urkm_RKM'*U;

%% Step 7)
clc;clear;
load data02.mat
K=3;
Rndstart=25;
Q=2;

%% Step 6)
[Urkm_RKMND,Arkm_RKMND, Yrkm_RKMND,frkm_RKMND,inrkm_RKMND] = REDKM(X, K, Q, Rndstart);
figure
gplotmatrix(Yrkm_RKMND,[],Urkm_RKMND*[1 2 3]','brkmcywv','+x^d.*os',[],'off','hist')