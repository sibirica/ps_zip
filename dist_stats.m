load('KarmaBY.mat', 'trat'); 
nt = size(trat,1);
st = 1;
R = zeros(nt,1);
ind = 1;
for t=1:nt
    frm = trat{t};
    R(t)=size(frm,1);
    for i=1:R(t)
        ds1 = zeros(R(t),1);
        ds2 = zeros(R(t),1);
        x1 = frm(i,1);
        y1 = frm(i,2);
        for j=1:R(t)
            x2 = frm(j,1);
            y2 = frm(j,2);
            if i==j
                ds1(j) = Inf;
            else
                ds1(j) = sqrt((x1-x2)^2+(y1-y2)^2);
            end
            if frm(i,3)+frm(j,3)==0 % opposite chiralities
                ds2(j)=ds1(j);
            else
                ds2(j)=Inf; % shouldn't be matched
            end
        end
        [dall(ind),~] = min(ds1);
        [ddiff(ind),~] = min(ds2);
        ind = ind+1;
    end
end
FS = 30;
TS = 26;
figure(2);
histogram(ddiff(st:end),20);
set(gcf,'color','w');
xlabel('$d$','FontSize',FS+2, 'Interpreter','latex')
ylabel('$N_d \ \times 10^{-3}$','FontSize',FS+2,'Interpreter','latex')
yticks(0:2000:8000)
xlim([0 145])
ylim([0 8500])
set(gca,'FontSize',TS+2,'FontName', 'Times')
pbaspect([3 2 1])
figure(1);
histogram(R(st:end),5:21);
set(gcf,'color','w');
xlabel('$n$','FontSize',FS, 'Interpreter','latex')
ylabel('$N_n$','FontSize',FS,'Interpreter','latex')
yticks(0:150:600)
xlim([5 21])
ylim([0 700])
set(gca,'FontSize',TS,'FontName', 'Times')
pbaspect([3 2 1])
disp("Mean PS separation:")
disp(mean(dall))
disp("Mean # PSs:")
disp(mean(R))
