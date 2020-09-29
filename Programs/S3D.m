figure
scatter(Y(:,1),Y(:,2),'filled')
xlabel('ACP Factor  1')
ylabel('ACP Factor 2')

figure
scatter(Ym50(:,1),Ym50(:,2))
xlabel('50 ObjCl Factor  1')
ylabel('50 ObjCl Factor 2')

figure
scatter(Ym30(:,1),Ym30(:,2))
xlabel('30 ObjCl Factor  1')
ylabel('30 ObjCl Factor 2')

figure
scatter(Ym20(:,1),Ym20(:,2))
xlabel('20 ObjCl Factor  1')
ylabel('20 ObjCl Factor 2')

%scatter3(Ym(:,1),Ym(:,2),Ym(:,3),'filled')

%zlabel('Factor 3')