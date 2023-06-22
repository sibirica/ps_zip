load('KarmaBY.mat', 'tracks');
ma = max(tracks(:,5)); % number of ids
lft = [];
displ = [];
trvl = [];
for id=1:ma
    trI = tracks(tracks(:,5)==id,:);
    lft = [lft size(trI,1)/53.33333];
    for j=2:size(trI,1)
        if trI(j,4)-trI(j-1,4)==1
           displ = [displ pdist([trI(j,1:2);trI(j-1,1:2)])];
        end
    end
    trvl = [trvl max(pdist(trI(:,1:2)))];
end
FS = 30;
TS = 26;
figure(1)
histogram(displ,20);
set(gca,'YScale','log')
set(gcf,'color','w');
xlim([0 8])
ylim([10 1e5])
xticks(0:2:8)
xlabel('{\it v}','FontSize',FS, 'Interpreter','latex')
ylabel('{\it N_v}','FontSize',FS, 'Interpreter','latex')
set(gca,'FontSize',TS,'FontName', 'Times')
pbaspect([3 2 1])
figure(2)
histogram(trvl,12);
set(gcf,'color','w');
xlabel('$r$','FontSize',FS, 'Interpreter','latex')
ylabel('$N_r$','FontSize',FS, 'Interpreter','latex')
xticks(0:40:120)
xlim([0 135])
set(gca,'FontSize',TS,'FontName', 'Times')
pbaspect([3 2 1])
set(gcf,'color','w');
figure(3)
histogram(lft,12);
set(gcf,'color','w');
xlabel('$l$','FontSize',FS, 'Interpreter','latex')
ylabel('$N_l$','FontSize',FS,'Interpreter','latex')
xticks(0:15:45)
xlim([0 45])
ylim([0 85])
set(gca,'FontSize',TS,'FontName', 'Times')
pbaspect([3 2 1])
set(gcf,'color','w');
disp('Median PS lifetime (in periods):')
median(lft)
disp('Mean PS velocity:')
mean(displ)
disp('Maximum PS displacement:')
max(trvl)