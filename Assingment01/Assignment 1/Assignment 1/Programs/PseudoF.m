function [Yo] = PseudoF(data)

% data: data matrix to clustering

eva  = evalclusters(data,'kmeans','Gap','ReferenceDistribution','PCA', 'KList',[2:10])
Yo   = eva.CriterionValues

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(Yo,'MarkerFaceColor',[0 0 1],'MarkerSize',5,'Marker','^',...
    'LineWidth',1,...
    'Color',[0 0 1]);

% Create xlabel
xlabel('Clusters');

% Create ylabel
ylabel('Pseudo-F function');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14,'XGrid','on','XTickLabel',...
    {'2','3','4','5','6','7','8','9','10'},'YGrid','on');
