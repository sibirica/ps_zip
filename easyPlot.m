function [] = easyPlot(data)
    mov = VideoWriter('original.avi','Motion JPEG AVI');
    mov.FrameRate = 24;
    open(mov);
    s = size(data);
    mini = min(data(:));
    maxi = max(data(:));
    colormap parula
    f = figure(1);
    set(gcf,'color','w');
    colorbar
    for t=1:s(3)
          imagesc(data(:,:,t), [mini maxi])
          axis square
          set(gca,'YDir','normal');
          axis off
          frame=getframe(f);
          writeVideo(mov,frame);
          if mod(t,15)==0
            cla; % speed 
          end
    end
    close(mov);
end