n = 4; %5
load('KarmaB.mat', 'nt');
dist = zeros(n,nt);
%total = dist;
actual = zeros(1,nt);
close all
figure
hold on
colors = 'kbrgcm';
for i=1:n
        filename = strcat('KarmaN1S',num2str(i));
        [~, dist(i,:), ~, actual(:)] = compare_track('KarmaB', filename, 0);
        nPS = unique(actual(1:end-1)); % last frame doesn't have PSs
        mdist = nPS;
        mind = 1;
        for j=nPS
            indj = (actual(:)==j);
            nums(mind) = sum(indj);
            mdist(mind) = sum(dist(i,indj))/(j*sum(indj));
            %if nums(mind)<20
            %    mdist(mind) = NaN;
            %end
            mind = mind+1;
        end
        s = strcat(colors(i),'o');
        d = strcat(colors(i),'-'); 
        plot(nPS, mdist, s);
        fitpoly = fit(nPS', mdist', 'poly2', 'Weight', nums);
        plot(fitpoly, d);
end