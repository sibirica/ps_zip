function [] = drawOnly(datapath, A, D1new, D2new, maskout, record, moviepath)
    % moviepath is the filename for saving the movie if record=1
    
    global nx nt W options
    
    mask3d = repmat(maskout,1,1,nt); % 3d version of mask
    A(mask3d==1) = NaN; % ignore data outside of ROI
    mini = min(min(min(A)));
    maxi = max(max(max(A)));
    
    % smoothed lvl plotting
    disp('Ready to plot!')
    warning('off','MATLAB:contour:ConstantData')
    colormap parula

    if record == 1
        mov = VideoWriter(moviepath,'Motion JPEG AVI');
        mov.FrameRate = 10;
        open(mov);
    end
    
    % if calculations already done, no need to repeat them; just draw.
    if exist(join([datapath 'X']),'file') == 2
        disp('Loading saved plotting data...')
        load(datapath, 'xs', 'ys', 'chi', 'numPS', 'charge'); % load necessary variables
        f = figure(1);
        for t=1:nt
           drawBare(A, mini, maxi, D1new, D2new, xs, ys, chi, t);
           txt = sprintf('# PS: %d; Charge: %d; t=%d',numPS(t),charge(t),t);
           title(txt)
           axis off
           drawnow;
           if record == 1
              frame=getframe(f);
              writeVideo(mov,frame);
           end
           if mod(t,15)==0
              cla; % speed 
           end
        end
        if record == 1
            close(mov);
        end
        return;
    end

    % distance from edge of mask determines what counts as ps
    em = edgemask(maskout); % create modified mask with edges
    [m1,m2] = find(em==1);
    maskdist = perform_fast_marching(W, [m1';m2'], options);
    dis = 4; % min allowed distance
    
    maxPS = 50; % maximum number of PS we could expect to get

    numPS = zeros(nt,1,'uint8');
    charge = numPS;
    xs = zeros(maxPS,nt,'single');
    ys = xs;
    chi = xs;
    C1 = zeros(2,8*nx,nt,'single');
    C2 = C1;
    f = figure(1);
    for t=1:nt
        % draw phase, lvl sets, and PSs (w/ chirality) and count them
        [num, ch, C1frame, C2frame, xi, yi, chiframe] = PSandDraw(A, mini, maxi, D1new, D2new, maskdist, t, dis);
        txt = sprintf('# PS: %d; Charge: %d; t=%d',num,ch,t);
        title(txt)
        numPS(t) = num;
        charge(t) = ch;
        C1(:,1:size(C1frame,2),t) = C1frame;
        C2(:,1:size(C2frame,2),t) = C2frame;
        xs(1:size(xi),t) = xi;
        ys(1:size(yi),t) = yi;
        chi(1:size(chiframe),t) = chiframe;
        axis off
        drawnow;
        % save for movie
        if record == 1
            frame=getframe(f);
            writeVideo(mov,frame);
        end
        if mod(t,20)==0
            cla % clearing the figure helps with speed
        end
    end
    save(datapath, 'numPS', 'charge', 'C1', 'C2', 'xs', 'ys', 'chi', '-append');
    done = 1;
    save(join([datapath 'X']), 'done'); % if this file exists the next time, PSandDraw will not run
    
    if record == 1
        close(mov);
    end
end

function mo2 = edgemask(maskout)
    mo2 = maskout;
    mo2(1,:) = 1;
    mo2(end,:) = 1;
    mo2(:, 1) = 1;
    mo2(:, end) = 1;
end

function [xc, yc] = contoursep(C)
    index = 1;
    xc = C(1,:);
    yc = C(2,:);
    while size(C,2)>index
        len = C(2,index); % number of points in this part
        xc(index) = NaN; % new separator
        yc(index) = NaN;
        index = index+len+1;
    end
end

function [num, ch, C1, C2, yi, xi, chi] = PSandDraw(A, mini, maxi, D1new, D2new, maskdist, t, dis)
    global nx
    num = 0;
    ch = 0;
    %imagesc(A(:,:,t), [mini maxi]);
    imagesc(A(:,:,t));
    hold on
    [C1, ~] = contour(D1new(:,:,t),[0,0],'-w','LineWidth',2);
    [C2, ~] = contour(D2new(:,:,t),[0,0],'-r','LineWidth',2);
    
    % separate contours into their component pieces
    [xc1, yc1] = contoursep(C1);
    [xc2, yc2] = contoursep(C2);
    
    % gradients of distance functions
    [gd1x, gd1y] = gradient(D1new(:,:,t));
    [gd2x, gd2y] = gradient(D2new(:,:,t));
    % find phase singularities
    if (size(xc1,2)>1 && size(xc2,2)>1)
        % for some reason x and y are switched!?
        [yi,xi] = intersections(xc1,yc1,xc2,yc2);
        % compute chirality for each singularity   
        ni = size(xi,1);
        chi = zeros(ni,1);
        for z=1:size(xi,1)
            % where is this intersection?
            xloc = xi(z);
            yloc = yi(z);
            neux = floor(xloc):ceil(xloc); % need meshgrid locally only
            neuy = floor(yloc):ceil(yloc);
            [Xnr,Ynr] = meshgrid(neux,neuy); % for local interpolation
            % skip this ps if it is close to the edge of mask
            %if interp2(Xnr,Ynr,maskdist(neux,neuy),xloc,yloc)<dis/nx
            if maskdist(round(xloc),round(yloc))<dis/nx
               continue
            end
            % find gradients
            gdh1x = interp2(Xnr,Ynr,gd1x(neux,neuy),xloc,yloc);
            gdh1y = interp2(Xnr,Ynr,gd1y(neux,neuy),xloc,yloc); 
            gdh2x = interp2(Xnr,Ynr,gd2x(neux,neuy),xloc,yloc);
            gdh2y = interp2(Xnr,Ynr,gd2y(neux,neuy),xloc,yloc); 
            chi(z) = gdh1x*gdh2y-gdh1y*gdh2x; % z-component of cross product
            if chi(z)>0 % CW rotation - positive chi - black
               % y and x switched back!?!?
               plot(yi(z),xi(z),'ko','MarkerSize',8,'MarkerFaceColor','k')
               num = num+1;
               ch = ch+1;
            elseif chi(z)<0 % CCW rotation - negative chi - white
               plot(yi(z),xi(z),'ko','MarkerSize',8,'MarkerFaceColor','w')
               num = num+1;
               ch = ch-1;
            else % we might get NaN next to the edge of the heart - grey
               plot(yi(z),xi(z),'ko','MarkerSize',8,'MarkerFaceColor','[0.5 0.5 0.5]')
            end 
        end
    end
end