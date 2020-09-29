function CalHar(data)

% data: data matrix to clustering

eva = evalclusters(data,'kmeans','CalinskiHarabasz','KList',[2:6])
Y   = eva.CriterionValues

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(Y,'MarkerFaceColor',[0 0 1],'MarkerSize',8,'Marker','^',...
    'LineWidth',2,...
    'Color',[0 0 1]);

% Create xlabel
xlabel('Clusters');

% Create ylabel
ylabel('Calinski-Harabasz function');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14,'XGrid','on','XTickLabel',...
    {'2','3','4','5','6'},'YGrid','on');
