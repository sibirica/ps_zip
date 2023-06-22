function []=sideBySide(filenameVm, filenameCa)
% Note: mini/maxi not included
    fileVm = sprintf('%s.mat', filenameVm);
    fileCa = sprintf('%s.mat', filenameCa);
    fileXVm = sprintf('%sX.mat', filenameVm);
    fileXCa = sprintf('%sX.mat', filenameCa);
    if (exist(fileXVm,'file') ~= 2) || ...
       (exist(fileXCa,'file') ~= 2)
        disp('You must finish processing the data first!')
        return
    end
    disp('Loading saved plotting data...')
    load(fileVm, 'A', 'D1new', 'D2new', 'nt', 'maskout');
    load(fileXVm,'xs', 'ys', 'xcs1', 'xcs2', 'ycs1', 'ycs2', 'chi', 'numPS', 'charge');
    AVm=A; D1Vm=D1new; D2Vm=D2new; x1Vm=xcs1; x2Vm=xcs2; y1Vm=ycs1; y2Vm=ycs2; xsVm=xs; ysVm=ys; chiVm=chi; numVm=numPS; chVm=charge;
    load(fileCa, 'A', 'D1new', 'D2new');
    load(fileXCa,'xs', 'ys', 'xcs1', 'xcs2', 'ycs1', 'ycs2', 'chi', 'numPS', 'charge');
    ACa=A; D1Ca=D1new; D2Ca=D2new; x1Ca=xcs1; x2Ca=xcs2; y1Ca=ycs1; y2Ca=ycs2; xsCa=xs; ysCa=ys; chiCa=chi; numCa=numPS; chCa=charge;
    mask3d = repmat(maskout,1,1,nt); % 3d version of mask
    AVm(mask3d==1) = NaN; % ignore data outside of ROI
    ACa(mask3d==1) = NaN;
    %miniVm = min(min(min(Vm)));
    %miniCa = min(min(min(Ca)));
    %maxiVm = max(max(max(Vm)));
    %maxiCa = max(max(max(Ca)));
    disp('Ready to plot!')
    
    clf;
    
    record = 1;
    if record == 1
        mov = VideoWriter('comparison.avi','Motion JPEG AVI');
        mov.FrameRate = 10;
        open(mov);
    end
    
    f = figure(1);
    colormap parula
    colorbar
    set(f, 'Position', [300 300 900 400]);
    for t=1:nt
        subplot(1,2,1);
        if mod(t,10)==0
            cla; % speed 
        end
        [C1Vm, C2Vm] = repackage(squeeze(x1Vm(:,t)),squeeze(y1Vm(:,t)),squeeze(D2Vm(:,:,t)));
        [C3Vm, C4Vm] = repackage(squeeze(x2Vm(:,t)),squeeze(y2Vm(:,t)),squeeze(D1Vm(:,:,t)));
        drawBare(AVm, xsVm, ysVm, chiVm, t,  C1Vm, '-w', C2Vm, ':w', C3Vm, ':k', C4Vm, '-k'); %miniVm, maxiVm
        txt = sprintf('method 1; # PS: %d; Charge: %d; t=%d',numVm(t),chVm(t),t);
        title(txt)
        subplot(1,2,2);
        if mod(t,10)==0
            cla; % speed 
        end
        [C1Ca, C2Ca] = repackage(squeeze(x1Ca(:,t)),squeeze(y1Ca(:,t)),squeeze(D2Ca(:,:,t)));
        [C3Ca, C4Ca] = repackage(squeeze(x2Ca(:,t)),squeeze(y2Ca(:,t)),squeeze(D1Ca(:,:,t)));
        drawBare(ACa, xsCa, ysCa, chiCa, t,  C1Ca, '-w', C2Ca, ':w', C3Ca, ':k', C4Ca, '-k'); %miniCa, maxiCa
        txt = sprintf('method 2; # PS: %d; Charge: %d',numCa(t),chCa(t));
        title(txt)
        if record == 1
             frame=getframe(f);
             writeVideo(mov,frame);
        end
    end
    if record == 1
        close(mov);
    end
end