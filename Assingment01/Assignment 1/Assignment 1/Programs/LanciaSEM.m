% Lancia SEM
pA=[1,1;2,1;3,1; 4,1];
%pC=[2,1;3,1]; % modello utilizzato suR e OpenMx
pC=[1,1;2,1;3,1];
pB=0;
pG=[1 1];
pSF=0;
pSZ=[1,1];
pSD=[1,1;2,2;3,3;4,4];
pSE=[1,1;2,2;3,3];
[AA,CC,GG,BB,SF,SZ,SD,SE,f]=SEM(X,Y,'S',pA,pC,pB,pG,pSF,pSZ,pSD,pSE,'OLS');
