function [] = drawBare(A, xs, ys, chi, t, mini, maxi, varargin)
    leastArg = 7; % least arguments accepted without any varargin
    numCurves = floor((nargin-leastArg)/2);
    imagesc(A(:,:,t), [mini maxi]);
    axis square
    set(gca,'YDir','normal');
    axis off
    set(gcf,'color','w');
    hold on
    for i=1:numCurves
        plotByPieces(varargin{2*i-1},varargin{2*i});
    end
    for z=1:size(xs,1)
        if chi(z,t)>0
            plot(xs(z,t),ys(z,t),'ko','MarkerSize',8,'MarkerFaceColor','k')
        elseif chi(z,t)<0
            plot(xs(z,t),ys(z,t),'ko','MarkerSize',8,'MarkerFaceColor','w')
        end
    end
    drawnow;
end