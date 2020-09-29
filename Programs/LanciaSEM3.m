% Lancia SEM3 solo analisi fattoriale confermativa
Y=[];
pA=[1,1;2,1;3,1;4,1; 5,2;6,2;7,2];
pC=0; %se pC=0 non ci sono variabili endogene non ci sono variabili latenti endogene
pB=0; %se non ci sono variabili endogene non ci sono relazioni fra variabili latenti endogene        
pG=[1,1;1,2];
pSF=[1,2]; %se SF=0 allora la matrice delle var lat esogene è identità 
pSZ=[1,1];
pSD=[1,1;2,2;3,3;4,4;5,5;6,6;7,7];
pSE=0;
[AA,CC,GG,BB,SF,SZ,SD,SE,f]=SEM(X,Y,'S',pA,pC,pB,pG,pSF,pSZ,pSD,pSE,'ML');
