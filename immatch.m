m = 4;
n = 5;
tscore = zeros(m,n);
dscore = zeros(m,n);
for i=1:m
    for j=1:n
        filename = strcat('LarmaN',num2str(i),'S',num2str(j));
        if exist(strcat(filename,'Y.mat'), 'file')
            [tscore(i,j), dscore(i,j)] = compare_new('KarmaB', filename);
        end
    end
end

% imagesc(tscore, [0 1])
% c = colorbar;
% set(gca,'ytick',1:4)
% set(gca,'xticklabel',{'256' '64' '32' '16' '8'})
% set(gca,'yticklabel',{'0' '0.1' '0.3' '1'})
% set(gcf,'color','w');
% figure(2)
% imagesc(dscore, [0, max(dscore(:))]);
% colorbar
% set(gca,'ytick',1:4)
% set(gca,'xticklabel',{'256' '64' '32' '16' '8'})
% set(gca,'yticklabel',{'0' '0.1' '0.3' '1'})
% set(gcf,'color','w');
tscore
dscore