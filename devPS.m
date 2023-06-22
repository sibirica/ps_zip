n = 5;
m = 4;
st = 500;
load('KarmaB.mat', 'nt');
%load('KarmaBX.mat', 'numPS');
load('KarmaBY.mat', 'trat');
nB = zeros(nt-st,1);
for i=st:nt-1  % beginning is unrepresentative, last frame doesn't have PSs
   nB(i-st+1) = size(trat{i},1); 
end
%nB = numPS(st:end-1);
nBu = unique(nB)'; 
RMSEsp = zeros(n-1,1);
close all
figure
hold on
colors = 'wkbrgcm'; % first one is currently ignored anyway
for i=2:n
        filename = strcat('LarmaN1S',num2str(i),'Y.mat'); %X.mat
        load(filename, 'trat'); %numPS
        numPS = zeros(nt-st,1);
        for k=st:nt-1
            numPS(k-st+1) = size(trat{k},1); 
        end
        RMSEsp(i-1) = rms(numPS-nB);
        stds = zeros(size(nBu));
        ind = 1;
        for j=nBu
            stds(ind) = std(numPS(nB(:)==j));
            ind = ind+1;
        end
        s = strcat(colors(i),'-o');
        plot(nBu, stds, s);
end
RMSEn = zeros(m-1,1);
figure
hold on
colors = 'wkbrgcm';
for i=2:m
        filename = strcat('LarmaN',num2str(i),'S1Y.mat'); %X.mat
        load(filename, 'trat'); %numPS
        numPS = zeros(nt-st,1);
        for k=st:nt-1
            numPS(k-st+1) = size(trat{k},1); 
        end
        RMSEn(i-1) = rms(numPS-nB);
        stds = zeros(size(nBu));
        ind = 1;
        for j=nBu
            stds(ind) = std(numPS(nB(:)==j));
            ind = ind+1;
        end
        s = strcat(colors(i),'-o');
        plot(nBu, stds, s);
end