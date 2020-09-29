% Classification
% Metodologia Statistica Avanzata 2012-2013
% Vichi

% Example iris data
load fisheriris
gscatter(meas(:,1), meas(:,2), species,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');
N = size(meas,1);

% linear and quadratic discriminat analysis
%                   sample      training set class_var               
ldaClass = classify(meas(:,1:2),meas(:,1:2),species);

% Compute the missclassification error
bad = ~strcmp(ldaClass,species);
ldaMISErr = sum(bad) / N


% Compute the confusiona matrix
[ldaResubCM,grpOrder] = confusionmat(species,ldaClass)


% compute the plot of the missclassified observations
pause
figure
gscatter(meas(:,1), meas(:,2), species,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');
hold on;
plot(meas(bad,1), meas(bad,2), 'kx');
hold off;


pause
figure
[x,y] = meshgrid(4:.1:8,2:.1:4.5);
x = x(:);
y = y(:);
j = classify([x y],meas(:,1:2),species);
gscatter(x,y,j,'grb','sod')

% Compute the quadratic discriminant regression
qdaClass = classify(meas(:,1:2),meas(:,1:2),species,'quadratic');

% Compute the missclassification error of the quadratic discriminant
% regression
bad = ~strcmp(qdaClass,species);
qdaMISErr = sum(bad) / N


pause
% Cross-validation
% first fix the random seed 

s = RandStream('mt19937ar','seed',0);
RandStream.setDefaultStream(s);


% use cvpartition to generate 10 disjoint stratified subsets.

cp = cvpartition(species,'k',10)


% Estimate the true test error for LDA using 10-fold stratified cross-validation.
ldaClassFun= @(xtrain,ytrain,xtest)(classify(xtest,xtrain,ytrain));
ldaCVMISErr  = crossval('mcr',meas(:,1:2),species,'predfun', ...
             ldaClassFun,'partition',cp)


% Estimate the true test error for QDA using 10-fold stratified cross-validation.
qdaClassFun = @(xtrain,ytrain,xtest)(classify(xtest,xtrain,ytrain,...
              'quadratic'));
qdaCVMISErr = crossval('mcr',meas(:,1:2),species,'predfun',...
           qdaClassFun,'partition',cp)

pause       

% The classregtree class creates a classification tree. 
t = classregtree(meas(:,1:2), species,'names',{'SL' 'SW' });



% It's interesting to see how the decision tree method divides the plane. 
% Use the same technique as above to visualize the regions assigned to each species.

[grpname,node] = t.eval([x y]);
gscatter(x,y,grpname,'grb','sod')
       


% draw a diagram of the decision rule and class assignments.
view(t); 



% Compute the missclassification error and the cross-validation error for classification tree.
dtclass = t.eval(meas(:,1:2));
bad = ~strcmp(dtclass,species);
dtMISErr = sum(bad) / N

dtClassFun = @(xtrain,ytrain,xtest)(eval(classregtree(xtrain,ytrain),xtest));
dtCVMISErr  = crossval('mcr',meas(:,1:2),species, ...
          'predfun', dtClassFun,'partition',cp)
      
% Pruning
% Try pruning the tree. First compute the MissErr for various of subsets of the original tree. 
% Then compute the cross-validation error for these sub-trees. 
% A graph shows that the MissErr is overly optimistic. It always decreases as the tree size grows, 
% but beyond a certain point, increasing the tree size increases the cross-validation error rate.
  
pause

MisErr = test(t,'resub');
[cost,secost,ntermnodes,bestlevel] = test(t,'cross',meas(:,1:2),species);
plot(ntermnodes,cost,'b-', ntermnodes,MisErr,'r--')
figure(gcf);
xlabel('Number of terminal nodes');
ylabel('Cost (misclassification error)')
legend('Cross-validation','MissClasErr')

[mincost,minloc] = min(cost);
cutoff = mincost + secost(minloc);
hold on
plot([0 20], [cutoff cutoff], 'k:')
plot(ntermnodes(bestlevel+1), cost(bestlevel+1), 'mo')

pt = prune(t,bestlevel);
view(pt)
legend('Cross-validation','MissClasErr','Min + 1 std. err.','Best choice')
hold off