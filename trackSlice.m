st = 260;
en = 315;
FS = 15;
TS = 11;
trackSec = tracks((tracks(:,4)>=st & tracks(:,4)<=en),:);
ids = unique(trackSec(:,5))';
for i=ids
   w = 1;
   [px, py, pt, ~] = path2(trackSec,i);
   if pt(end)<en || pt(1)>st
      w = 2;
   end
   plot3(px, py, pt, 'LineWidth', w);
   set(gca,'xtick',[1 64:64:256])
   set(gca,'ytick',[1 64:64:256])
   set(gca,'ztick',[270 290 310])
   set(gca,'FontSize',TS,'FontName', 'Times')
   xlabel('$x$','FontSize',FS,'Interpreter','latex')
   ylabel('$y$','FontSize',FS,'Interpreter','latex')
   zlabel('$t$','FontSize',FS,'Interpreter','latex')
   axis([1 256 1 256 st en])
   hold on
end
set(gcf,'color','w');