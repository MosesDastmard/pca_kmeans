% Linear discriminant analysis Example
load fisheriris

% 150 x 2 
% Three classes 'setosa', 'versicolor', 'virginica'

gscatter(meas(:,1), meas(:,2), species,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');
N = size(meas,1);
uk=ones(3,1);


% linear discriminat analysis

ldaClass = classify(meas(:,1:2),meas(:,1:2),species);

% confusion matrix

[Colda,grpOrder] = confusionmat(species,ldaClass)

pause 

ClErLDA = (uk'*Colda*uk - uk'*diag(Colda))/N

hold on;
bad = ~strcmp(ldaClass,species);
plot(meas(bad,1), meas(bad,2), 'kx');
hold off;

pause

[x,y] = meshgrid(4:.1:8,2:.1:4.5);
x = x(:);
y = y(:);
j = classify([x y],meas(:,1:2),species);
gscatter(x,y,j,'grb','sod')
   
pause

% Cross validation k-folds

s = RandStream('mt19937ar','seed',0);
RandStream.setDefaultStream(s);

cp = cvpartition(species,'k',10)

pause

ldaClassFun= @(xtrain,ytrain,xtest)(classify(xtest,xtrain,ytrain));
ldaCVErr  = crossval('mcr',meas(:,1:2),species,'predfun', ...
             ldaClassFun,'partition',cp)

         
pause
% quadratic discriminant analysis


qdaClass = classify(meas(:,1:2),meas(:,1:2),species,'quadratic')

[Coqda,grpOrder] = confusionmat(species,qdaClass)

pause 

ClErQDA = (uk'*Coqda*uk - uk'*diag(Coqda))/N

hold on;
bad = ~strcmp(qdaClass,species);
plot(meas(bad,1), meas(bad,2), 'kx');
hold off;

pause

% Cross validation k-folds

s = RandStream('mt19937ar','seed',0);
RandStream.setDefaultStream(s);

cp = cvpartition(species,'k',10)
% 
qdaClassFun = @(xtrain,ytrain,xtest)(classify(xtest,xtrain,ytrain,...
              'quadratic'));
qdaCVErr = crossval('mcr',meas(:,1:2),species,'predfun',...
           qdaClassFun,'partition',cp)