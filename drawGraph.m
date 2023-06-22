function [] = drawGraph(grps, tracks, xs, ys, chi, t, maskout)
    %style = {'-w', ':w', ':k', '-k', '-r', ':r', ':y', '-y'};
    style = {'-w', ':w', ':k', '-k'};
    figure(6)
    clf
    imagesc(maskout)
    axis square
    set(gca,'YDir','normal');
    axis off
    set(gcf,'color','w');
    hold on
    grp = grps{t};
    tr = tracks(tracks(:,3)==t,:);
    for z=1:size(xs,1)
        if chi(z,t)>0
            plot(xs(z,t),ys(z,t),'ko','MarkerSize',8,'MarkerFaceColor','k')
        elseif chi(z,t)<0
            plot(xs(z,t),ys(z,t),'ko','MarkerSize',8,'MarkerFaceColor','w')
        end
    end
    for i=1:size(grp,1)
        type = grp(i,1); id1 = grp(i,2); id2 = grp(i,3);
        q1 = tr(tr(:,4)==id1,:);
        q2 = tr(tr(:,4)==id2,:);
        plot([q1(1)+(type-2)/4 q2(1)],[q1(2) q2(2)+(mod(type,2)-0.5)],style{type},'LineWidth',2);
    end
end