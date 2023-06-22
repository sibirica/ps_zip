function [x,y] = pick_up_a_trace_NEW(DATA, frame)
if nargin == 1
    frame = 20;
end
figure(99),
subplot(232);
imagesc(mat2gray(DATA(:,:,frame))),colormap('gray'), axis square
button =0; drawnow
while button ~=3
    subplot(232);
    drawnow
    [x,y,button]=ginput(1);
    x=round(x); y=round(y);

    if button == 32 
        I(x,y) = max(max(DATA(:, :, frame)));
        imagesc(I);
    end

    X = squeeze(mean(mean(DATA(y-1:y+1, x-1:x+1, :))));
    subplot(2,3,4:6)
    plot(X);
    grid on; grid minor;
    title(strcat('row: ', num2str(y), ' col: ', num2str(x)))
    drawnow
    pause(0.05);
end
close(99);

