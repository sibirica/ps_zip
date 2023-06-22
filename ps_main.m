function [fail] = ps_main(opts, filename, datain, maskfile)
% from opt struct:
% exp is 1 for optical mapping data, 2 for data extracted from movie,
% 3 for data recorded on a sphere (i.e. basket catheter)
% fit is 1 if time derivatives calculated by least squares (much slower)
% record is 1 if saving a movie
% trQ is 1 if tracking is repeated
% mode is 1 for level sets defined by phase, 2 for level sets
% defined by max/min/back/front, 3 for phase + sparsity (no level sets)
% sc is the degree to which the data is spatially sparsified
% noise - percentage noise added artificially to A_hipa
% important parameters in this method can also be set
% (old option) fou is 1 for Fourier filtering of distance function
% filename is where calculations will be saved afterwards
% (optional) datain is the raw data; if provided, no calculations skipped
% (optional) maskfile is the mat file containing the mask
% 1) Conduction velocity
% 2) Front curvature
% 3) Detect PSs exiting the domain
% 4) Add support for spherical configuration: fix interpolation & QRS filter
% 5) 3D version of algorithm
    global hp dp tp
    hp = 1.57079633; % half pi
    dp = 4.71238898; % three pi / 2
    tp = 6.28318530; % two pi
    
    
    global record trQ moviepath hq
    if isfield(opts,'mode')
        mode = opts.mode;
    else
        mode = 2;
    end
    if isfield(opts,'exp')
        exp = opts.exp;
    else
        exp = 1;
    end
    if isfield(opts,'record')
        record = opts.record;
    else
        record = 1;
    end
    if isfield(opts,'trQ')
        trQ = opts.trQ;
    else
        trQ = 1;
    end
    if isfield(opts,'hq')
        hq = opts.hq;
    else
        hq = 0;
    end
    if isfield(opts,'moviepath')
        moviepath = opts.moviepath;
    else
        moviepath = filename;
    end
    
    addpath(genpath(pwd))
    file = sprintf('%s.mat', filename);
    fileX = sprintf('%sX.mat', filename);
    % load data if not provided
    if nargin<3
        if exist(fileX, 'file') == 2
            fail = TopandPlot(filename);
            return
        elseif exist(file, 'file') == 2 % calculations already saved
            load(file, 'mode');
            if mode==3
                findPSsparse(filename);
                fail = TopandPlot(filename);
            else
                load(file, 'D1new', 'D2new', 'maskout', 'nx', 'ny', 'nt'); % this is all we might need
                fail = PStrack(filename, D1new, D2new, maskout);
            end
            return % that is all that needs to be done, so we can exit
        end
        data = single(opread('E:\PS\raw data\pig\Rec003.dat'));
        %load('C:\Users\Daniel\Dropbox\PS\daniel_yay\Vm.mat')
        %data = Vm;
    else
        data = single(datain);
    end
    
    global nx ny nt
    [nx,ny,nt] = size(data);

    if isfield(opts,'fit')
        fit = opts.fit;
    else
        fit = 0;
    end
    if isfield(opts,'sc')
        sc = opts.sc;
    else
        sc = 0;
    end
    if isfield(opts,'noise')
        noise = opts.noise;
    else
        noise = 0;
    end
    if isfield(opts,'MPP') % minimum peak prominence (for mode 1 maxima)
       MPP = opts.MPP;
    else
       MPP = 0.05;
    end
    if isfield(opts,'MPP1') % for peaks of signal (mode 2)
       MPP1 = opts.MPP1;
    else
       if noise==0 && opts.exp~=1 
            MPP1 = 0.01;
       elseif noise<=0.1
            MPP1 = 0.05;
       else
            MPP1 = 0.25; % high noise
       end
    end
    if isfield(opts,'MPP2') % for peaks of derivative (mode 2)
       MPP2 = opts.MPP2;
    else
       if noise==0 && opts.exp~=1 
            MPP2 = 0.05;
       elseif noise<=0.1
            MPP2 = 0.1;
       else
            MPP2 = 0.25; % high noise
       end
    end
    if isfield(opts,'sig')
        sig = opts.sig;
    elseif exp==1
        sig = [2 2 3]; % gaussian width for smoothing
    elseif exp==2 % JPEG MODE / artificial noise
        %sig = [1.5 1.5 3];
        sig = [1e-6 1e-6 5];
    elseif exp==3
        sig = [1e-6 1e-6 10];
    end 
    if isfield(opts,'sig2')
       sig2 = opts.sig2;
    else
       sig2 = [1e-6 1e-6 1e-6]; % [1 1 2], [1.5 1.5 2+(sc>0)], [1e-6 1e-6 3]
    end
    if isfield(opts,'sigD')
       sigD = opts.sigD;
    elseif exp==1
        %sigD = 4;
        sigD = max(nx,ny)/32;
    else
        %sigD = 2;
        sigD = max(nx,ny)/64;
    end

    % preprocessing
    disp('Data preprocessing...')
    % for Rec003/012 only (fix data glitches):
    if exp==1
        %optical mapping
        %data(1,1,:) = data(2,1,:);
        %data(128,:,:) = data(127,:,:);
        data = -data;
    end
    st = 1; % first frame of recording processed: 003:5,012:200,other:1
    if mod(nt-st+1,2)==1
       st = st+1; 
    end
    L = nt-st+1; % # of frames analyzed; can be set manually
    %L = 2100; % 003:4996,012:2100
    if exp==1
        de = ceil(L/500); % # of lowest frequencies thrown out
    % may already be processed earlier
    %elseif exp==3
    %   de = ceil(L/2000);
    else
        de = 1;
    end
    filter = (exp==1); % || exp==3); % see above
    [A_hipa, ffa, lp] = procRawData(data,st,L,de,filter);
    T = median(squeeze(reshape(lp,[nx*ny,1]))); % estimated period is median of local periods
    nt = size(A_hipa,3);
    
    % compute mask
    if nargin>3
        load(sprintf('%s.mat', maskfile)) % load existing mask
        thres = NaN;
    elseif exp==1
        [maskout, thres] = makeMask(ffa); 
    else
        thres = NaN;
        maskout = zeros(nx,ny);
    end
    % constants
    if exp==1 || exp==3
        deltat = uint16(T/5);
        cut = max(30,de); % # of discarded frames at beginning and end (Gibbs)
        A_hipa = A_hipa(:,:,1+cut:end-cut);
        nt = nt-2*cut;
    else
        deltat = uint16(T/5);
    end

    % normalization per pixel
    if exp==1 || exp==3
        for x=1:nx
            for y=1:ny
                mi = min(A_hipa(x,y,:));
                ma = max(A_hipa(x,y,:));
                if ma-mi==0
                   maskout(x,y)=1;
                   A_hipa(x,y,:) = 0;
                else
                   A_hipa(x,y,:) = (A_hipa(x,y,:)-mi)/(ma-mi);
                end
            end
        end
    else % or globally
        mi = min(A_hipa(:));
        ma = max(A_hipa(:));
        A_hipa(:) = (A_hipa(:)-mi)/(ma-mi);
    end
    if noise>0
        A_hipa = A_hipa+noise*randn(size(A_hipa));
    end
        
    % sparsify
    if sc>0
        disp('Sparsifying...')
        [A_hipa, As] = sparsify(A_hipa,sc);
    elseif mode==3
        As = A_hipa;
    end
    
    % smooth data
    disp('Smoothing data...')
    if exp>0
        A = imgaussfilt3(A_hipa,sig);
    else
        A = A_hipa;
    end
    mask3d = repmat(maskout,1,1,nt); % 3d version of mask
    %A(mask3d==1) = NaN;

    if mode==1
        % locate maxima (phase = 0)
        disp('Locating maxima...')
        maxiind = findmaxiind(A,maskout,deltat,MPP);

        % initialize phase at all points and times
        disp('Initializing phase...')
        ph = makePhase(maxiind,maskout);

        % detect phase crossing 0, pi/2, pi, 3pi/2
        disp('Building level sets...')
        lvl = makelvl(ph,maskout);
    elseif mode==3
        c = 1.7; % TO DO: find automatically
        omega = 2*pi/T;
        k = omega/c;
        disp('Locating maxima...')
        maxiind = findmaxiind(As,maskout,deltat,MPP,sc);

        % initialize phase at all points and times
        disp('Initializing phase...')
        ph = makePhase(maxiind,maskout,sc);
    else   
        % find Adiff (for wavefront/back), Adiff2 (classifying signed
        % distance)
        disp('Calculating time derivatives of data...');
        if fit==1
            twind = 11;
            Adiff = zeros(size(A));
            Adiff2 = Adiff;
            w = waitbar(0,'Differentiating...'); 
            for x=1:nx
                waitbar(x/nx);
                for y=1:ny
                    [Adiff(x,y,:), Adiff2(x,y,:)] = diff2(squeeze(A(x,y,:)), twind);
                end
            end
            close(w);
            Adiff = Adiff/max(Adiff(:));
            Adiff2 = Adiff2/max(Adiff2(:));
        else
            %Adiff = diff(A,1,3);
            %with central differencing (not sure it's better though)
            Adiff = tdiff(A);
            if exp~=0
                Adiff = imgaussfilt3(Adiff, sig2);
            end
            Adiff = Adiff/max(Adiff(:)); % useful to normalize
            %Adiff(:,:,nt)=Adiff(:,:,nt-1); 
            %Adiff2 = diff(Adiff,1,3);
            Adiff2 = tdiff(Adiff);
            if exp~=0
                Adiff2 = imgaussfilt3(Adiff2, sig2);
            end
            Adiff2 = Adiff2/max(Adiff2(:));
            %Adiff2 = cat(3, Adiff2(:,:,1), Adiff2);
        end
        disp('Finding peaks...');
        maxi1 = findmaxiind(A, maskout, deltat, MPP1);
        maxi2 = findmaxiind(-A, maskout, deltat, MPP1);
        maxi3 = findmaxiind(Adiff, maskout, deltat, MPP2);
        maxi4 = findmaxiind(-Adiff, maskout, deltat, MPP2);
        estT = estperiod(maxi1);
        disp('Regularizing derivatives...');
        Adiff = regularize(Adiff, maskout, maxi1, maxi2);
        Adiff2 = regularize(Adiff2, maskout, maxi3, maxi4);
        disp('Building level sets...');
        lvl = makelvl2(Adiff,Adiff2,maskout);
    end
    if mode~=3
        % fast marching
        disp('Computing distance functions...')
        [D1, D2] = dist(lvl);
        disp('Making distance functions signed...')
        if mode==1
            D1s = D1.*sign(sin(ph));
            D2s = D2.*sign(cos(ph));
        else
            D1s = -D1.*sign(Adiff);
            D2s = -D2.*sign(Adiff2);
        end
        % smooth distance function by removing high spatial frequencies
        disp('Smoothing distance functions...')
        D1new = imgaussfilt(D1s, sigD);
        D2new = imgaussfilt(D2s, sigD);
        D1new(mask3d==1) = NaN; % don't need this anymore outside of ROI
        D2new(mask3d==1) = NaN;
        % time derivative of distance function 1 (smoothed)
        dD = imgaussfilt(tdiff(D1new), sigD);
        cv = abs(dD); % conduction velocity %./estT
    end
    disp('Saving data...')
    varN = {'mode', 'exp', 'sc', 'noise', 'maskout', 'data', 'A', 'nx', 'ny', 'nt', ...
        'T', 'lvl'};
    if mode~=3
        varN = horzcat(varN, {'D1new', 'D2new'});
    end
    if mode==1
        varN = horzcat(varN, {'MPP', 'ph'});
    elseif mode==2
        varN = horzcat(varN, {'MPP1', 'MPP2', 'estT', 'dD', 'cv'}); 
        % , 'Adiff', 'maxi1', 'maxi2', 'maxi3', 'maxi4'
    else
        varN = horzcat(varN, {'ph', 'As', 'k'});
    end
    if exp==1
        varN = horzcat(varN, {'lp', 'thres'});
    end
    if exp>=1
        varN = horzcat(varN, {'A_hipa', 'sig'});
    end
    save(file, varN{:}); % save named variables
    % call function for next part
    if mode==3
        findPSsparse(filename);
        fail = TopandPlot(filename);
    else
        fail = PStrack(filename, D1new, D2new, maskout);
    end
end

function [fail] = PStrack(filename, D1new, D2new, maskout)
    global nt
    
    fileX = sprintf('%sX.mat', filename);

    % distance from edge of mask determines what counts as ps
    em = edgemask(maskout); % create modified mask with edges
    maskdist = bwdist(em);
    dis = 2; % min allowed distance
    numPS = zeros(nt,1);
    charge = numPS;
    xs = [];
    ys = xs;
    Is = xs; % segment on curve1 for intersection point
    Js = xs;
    chi = xs;
    xcs1 = [];
    ycs1 = xcs1;
    xcs2 = xcs1;
    ycs2 = xcs1;
    f = figure(1);
    disp('Finding phase singularities...');
    w = waitbar(0,'Finding phase singularities...');
    for t=1:nt
        if mod(t,50)==0
           waitbar(t/nt); 
        end
        [num, ch, xc1, yc1, xc2, yc2, xi, yi, chiframe, I, J] = findPS(D1new, D2new, maskdist, t, dis);
        numPS(t) = num;
        charge(t) = ch;
        s1 = size(xc1,2);
        s2 = size(xc2,2);
        xcs1(1:s1,t) = xc1;
        ycs1(1:s1,t) = yc1;
        xcs2(1:s2,t) = xc2;
        ycs2(1:s2,t) = yc2;
        xs(1:size(xi),t) = xi;
        ys(1:size(yi),t) = yi;
        Is(1:size(I),t) = I;
        Js(1:size(J),t) = J;
        chi(1:size(chiframe),t) = chiframe;
    end
    close(w);
    disp('Saving data...')
    % if this file exists the next time, findPS will not run
    save(fileX, 'numPS', 'charge', 'xcs1', ... 
    'ycs1', 'xcs2', 'ycs2', 'xs', 'ys', 'Is', 'Js', 'chi');

    fail = TopandPlot(filename);
end

function [fail] = TopandPlot(filename)
    bg = 0; % use A/A_hipa or data

    file = sprintf('%s.mat',filename);
    fileX = sprintf('%sX.mat',filename);
    fileY = sprintf('%sY.mat',filename);
    disp('Loading saved plotting data...')
    global record trQ moviepath hq nx ny T
    load(file, 'exp', 'mode', 'maskout', 'nx', 'ny', 'nt', 'T');
    if mode~=3
        load(file, 'D1new', 'D2new');
    end
    if mode==1
        load(file, 'ph');
    elseif ~bg
        %if exp>=1
        %   load(file, 'A_hipa');
        %else
           load(file, 'A');
        %end
    else
        load(file, 'data');
        %load(file, 'cv');
    end
    load(fileX, 'xs', 'ys', 'chi', 'numPS', 'charge'); 
    if mode~=3
        load(fileX, 'xcs1', 'xcs2', 'ycs1', 'ycs2', 'Is', 'Js');
    end
    if trQ==1
        sv = 14; %20 - most stuff, or 14
        disp('Tracking PSs!')
        [tracks, L, vmax] = trackPS(xs,ys,chi,sv);
    else
        load(fileY, 'tracks', 'L', 'vmax')
    end
    if numel(tracks)==0
       fail = 1;
    else
       fail = 0;
    end
    mask3d = repmat(maskout,1,1,nt); % 3d version of mask
    if mode==1
        A_hipa = ph;
        A_hipa(mask3d==1) = NaN;
    elseif ~bg
        %if exp<1
            A_hipa = A;
        %end
        A_hipa(mask3d==1) = NaN; % ignore data outside of ROI
    else
        A_hipa = data;
        %A_hipa = cv;
    end
    hack = 0;
    mini = min(A_hipa(:))-hack; % colorbar hack
    maxi = max(A_hipa(:))+hack;
    %maxi = median(A_hipa(:))*10;
    %mini = 0; maxi = 1;

    if record == 1
        quality = {'Motion JPEG AVI', 'Uncompressed AVI'};
        mov = VideoWriter(moviepath,quality{hq+1});
        mov.FrameRate = 24;
        open(mov);
        disp('Ready to plot!')
        warning('off','MATLAB:contour:ConstantData')
        clf;
        colormap parula
        f = figure(1);
        %rng = (maxi-mini)/5;
        %locs = round(mini,1):round(rng,1):round(maxi,1);
        %lbls = {};
        %for i=1:length(locs)
        %    lbls{i} = sprintf('%0.1f',locs(i));
        %end
        %colorbar('northoutside','Ticks', locs, 'TickLabels', lbls,'FontSize',20)
        colorbar % often doesn't show up for some reason!
    end

    if mode~=3
        grps = cell(nt,1);
    end
    trat = cell(nt,1);
    for t=1:nt
       if ~fail
           trat{t} = tracks(tracks(:,4)==t,:);
       end
       if mode~=3
           [C1, C2, grp1] = repackage(squeeze(xcs1(:,t)),squeeze(ycs1(:,t)), ...
           squeeze(D2new(:,:,t)),trat{t},xs(:,t),ys(:,t),Is(:,t),0,mode);
           [C3, C4, grp2] = repackage(squeeze(xcs2(:,t)),squeeze(ycs2(:,t)), ...
           squeeze(D1new(:,:,t)),trat{t},xs(:,t),ys(:,t),Js(:,t),1,mode);
           grp = vertcat(grp1,grp2);
           grps{t} = grp;
       end
       if record==1
          if mode==3
              drawBare(A_hipa, xs, ys, chi, t, mini, maxi);
          else
              drawBare(A_hipa, xs, ys, chi, t, mini, maxi, C1, '-w', C2, ':w', C3, ':k', C4, '-k');
          end
          txt = sprintf('# PS: %d; Charge: %d; t=%d',numPS(t),charge(t),t);
          title(txt)
          frame=getframe(f);
          writeVideo(mov,frame);
          if mod(t,15)==0
            cla; % speed 
          end
       %elseif mod(t,100)==0
       %    t
       end
    end
    if record == 1
        close(mov);
    end
    if ~fail
        if mode==3
            disp('Saving...')
            save(fileY,'tracks','L','vmax','trat');
        else
            disp('Analyzing transitions...')
            [cre, ane, typeCr, typeAn] = ps_tran(tracks, grps, vmax, nt);
            cre = horzcat(typeCr, cre);
            ane = horzcat(typeAn, ane);
            disp('Saving...')
            save(fileY, 'tracks', 'L', 'vmax', 'trat', 'grps', 'cre', 'ane');
        end
    end
end

function [A_hipa, ffa, lp] = procRawData(data,st,L,de,filter)
    global nx ny
    ce = L/2; % center
    A_hipa = zeros(nx,ny,L,'single');
    ffa = zeros(nx,ny,ce,'single');
    w = waitbar(0,'Data preprocessing...');
    for x=1:nx
        waitbar(x/nx);
        for y=1:ny
            Vmx = squeeze(data(x,y,st:st+L-1));
            ff = fft(Vmx); % fft of signal at a point
            % find local period
            ffx = ff(1:ce); % consider first half only
            ffx(1:de) = 0; % ignore low frequencies
            ffa(x,y,:) = abs(ffx);

            if filter % apply high-pass filter to data
                ffs = fftshift(ff);
                ffs(ce-de:ce+de) = 0;
                A_hipa(x,y,:) = real(ifft(ifftshift(ffs)));
            end
        end
    end
    [~, lc] = max(ffa,[],3); % max along time dimension
    lp = L./lc;
    if ~filter
        A_hipa = data(:,:,st:st+L-1);
    end
    close(w);
end

function [maskout, thres] = makeMask(ffa) %lp,T
    % Procedure: first, click a few points inside the heart close to the
    % heart boundary and spread along it, and press enter. 
    % Then, click a few outside, and enter again.
    % If that looks impossible, immediately press enter and 
    % an alternative mask will be used.
    global nx ny
    disp('We need your help to make a mask!')
    %V = var(A_hipa,0,3);
    %mV = max(V(:));
    specmeans = zeros(nx,ny);
    for x=1:nx
        for y=1:ny
            specmeans(x,y) = mean(squeeze(ffa(x,y,:)));
        end
    end
    mm = max(specmeans(:));
    fm = figure(12);
    it = specmeans/mm;
    imagesc(it);
    set(gca,'YDir','normal');
    disp('Click a few times inside the heart and press enter.')
    [Xin,Yin]=ginput;
    nin = size(Xin);
    if nin==0
       %maskoutR = (lp<T/2 | lp>T*2);
       %sms = 0.8; % smoothing sigma (this approach is very broken!!!)
       %maskout = logical(imgaussfilt(double(maskoutR),sms));
       maskin = roipoly(it);
       set(gca,'YDir','normal');
       %imagesc(maskoutR);
       maskout = ~maskin;
       thres = NaN;
       close(fm);
       return
    end
    disp('OK! Click a few times outside and press enter again.')
    [Xout,Yout] = ginput;
    nout = size(Xout);
    guessIn = zeros(nin);
    for i=1:nin
       guessIn(i) = specmeans(round(Xin(i)),round(Yin(i)))/mm; 
    end
    minIn = min(guessIn);
    guessOut = zeros(nout);
    for i=1:nout
       guessOut(i) = specmeans(round(Xout(i)),round(Yout(i)))/mm; 
    end
    maxOut = max(guessOut);
    thres = (minIn+maxOut)/2;
    maskout = (specmeans/mm<thres);
    close(12);
end

function maxiind = findmaxiind(A,maskout,deltat,MPP,sc)
    global nx ny nt
    w = waitbar(0,'Locating maxima...');
    if nargin<5
        nxs = nx;
        nys = ny;
    else
        nxs = nx/2^sc+(sc>0);
        nys = ny/2^sc+(sc>0);
    end
    maxiind = zeros(nxs,nys,uint16(nt/deltat));
    for x=1:nxs
        waitbar(x/nxs);
        for y=1:nys
            if maskout(x,y)
                continue
            end
            [~, maxs] = findpeaks(squeeze(A(x,y,:)),'MinPeakDistance',deltat,'MinPeakProminence',MPP);
            maxiind(x,y,1:size(maxs,1)) = maxs;
        end
    end
    close(w);
end

function estT = estperiod(maxiind)
    global nx ny nt
    periods = diff(maxiind, [], 3);
    estT = zeros(nx, ny, nt);
    for x=1:nx
        for y=1:ny
            maxs = maxiind(x,y,:);
            nMaxs = nnz(maxs);
            estT(x,y,1:maxs(1)) = periods(x,y,1);
            for i=1:nMaxs-1
                estT(x,y,maxs(i):maxs(i+1)) = periods(x,y,i);
            end
            if nMaxs>1
                estT(x,y,maxs(nMaxs):end) = periods(x,y,nMaxs-1);
            else
                estT(x,y,:) = NaN;
            end
        end
    end
end

function ph = makePhase(maxiind, maskout, sc)
    % phase varies linearly from 0 to 2*pi (tp) between maxima
    global nx ny nt tp
    if nargin<3
        nxs = nx;
        nys = ny;
    else
        nxs = nx/2^sc+(sc>0);
        nys = ny/2^sc+(sc>0);
    end
    ph = zeros(nxs,nys,nt,'single'); % phase
    w = waitbar(0,'Computing phase...');
    for x=1:nxs
        waitbar(x/nxs);
        for y=1:nys
            if maskout(x,y)
                continue
            end
            ind = 1;
            indi = maxiind(x,y,ind);
            indf = maxiind(x,y,ind+1);
            diff = indf-indi;
            oldindi = indi-diff;
            for t=1:indi-1
                p = mod((t-oldindi)/diff,1);
                ph(x,y,t) = tp*p; % phase before max #1
            end
            while(indf~=0)
                for t=indi:indf-1
                   p = (t-indi)/diff;
                   ph(x,y,t) = tp*p;
                end
                ind = ind+1;
                indi = maxiind(x,y,ind);
                indf = maxiind(x,y,ind+1);
                diff = indf-indi;
            end
            if(ind==1)
                continue
            end
            indf = indi;
            indi = maxiind(x,y,ind-1);
            diff = indf-indi;
            for t=indf:nt
                p = mod((t-indf)/diff,1);
                ph(x,y,t) = tp*p;
            end
        end
    end
    close(w);
end

function Ar = regularize(A, maskout, maxs, mins)
    global nx ny
    for x=1:nx
        for y=1:ny
            if maskout(x,y)
                continue
            end
            la = find(maxs(x,y,:)==0,1)-1;
            li = find(mins(x,y,:)==0,1)-1;
            if isempty(la)
                la = size(maxs,3);
            end
            if isempty(li)
                li = size(mins,3);
            end
            st = 1;
            inda = 1;
            indi = 1;
            inc0 = (maxs(x,y,1)<mins(x,y,1));
            sgn = 2*(inc0-0.5);
            while inda<=la && indi<=li
                ma = maxs(x,y,inda);
                mi = mins(x,y,indi);
                if inda<la
                    na = maxs(x,y,inda+1);
                else
                    na = Inf;
                end
                if indi<li
                    ni = mins(x,y,indi+1);
                else
                    ni = Inf;
                end
                inc = (ma<mi);
                if inc>inc0
                    indi = indi-1;
                    if inda>1
                        mi = floor(0.5*(maxs(x,y,inda-1)+ma));
                    else
                        mi = floor(0.5*ma);
                    end
                elseif inc0>inc
                    inda = inda-1;
                    if indi>1
                        ma = floor(0.5*(maxs(x,y,indi-1)+mi));
                    else
                        ma = floor(0.5*mi);
                    end
                elseif na < mi
                    indi = indi-1;
                    mi = floor(0.5*(ma+na));
                elseif ni < ma
                    inda = inda-1;
                    ma = floor(0.5*(mi+ni));
                end
                mid = min(ma,mi);
                en = max(ma,mi);
                % this is a horrible hack to prevent sign(A) = 0
                A(x,y,st:mid) = (eps+abs(A(x,y,st:mid)))*sgn;
                A(x,y,mid+1:en) = -(eps+abs(A(x,y,mid+1:en)))*sgn;
                inda = inda+1;
                indi = indi+1;
                st = en+1;
            end
            if inda-la+indi-li==1 && la>0 && li>0 % one extra point
               mae = maxs(x,y,la);
               mie = mins(x,y,li);
               ince = (mae<mie);
               mid = max(mae, mie);
               if inc0 == ince % edge case such as ... mi ma ma
                  extra = floor((en+mid)/2);
                  A(x,y,st:extra) = (eps+abs(A(x,y,st:extra)))*sgn;
                  A(x,y,extra+1:mid) = -(eps+abs(A(x,y,extra+1:mid)))*sgn;
                  A(x,y,mid+1:end) = (eps+abs(A(x,y,mid+1:end)))*sgn;
               end
               A(x,y,st:mid) = (eps+abs(A(x,y,st:mid)))*sgn;
               A(x,y,mid+1:end) = -(eps+abs(A(x,y,mid+1:end)))*sgn;
            else
               A(x,y,st:end) = (eps+abs(A(x,y,st:end)))*sgn;
            end
        end
    end
    Ar = A;
end

function lvl = makelvl(ph, maskout)
    global nx ny nt hp dp tp
    lvl = zeros(nx,ny,nt,'single'); % level sets
    w = waitbar(0,'Building level sets...');
    for x=2:nx
        waitbar(x/(nx-1));
        for y=2:ny
            if maskout(x,y)
                continue
            end
            for t=1:nt
                a = ph(x-1,y,t);
                b = ph(x,y-1,t);
                H = ph(x,y,t);
                aa = mod(a+pi,tp)-pi; % for 0, tp problem
                bb = mod(b+pi,tp)-pi;
                HH = mod(H+pi,tp)-pi;
                if(aa*HH<=0 || (a-pi)*(H-pi)<0)
                    lvl(x,y,t) = 1;
                    lvl(x-1,y,t) = 1;
                elseif((a-hp)*(H-hp)<0 || (a-dp)*(H-dp)<0)
                    lvl(x,y,t) = -1;
                    lvl(x-1,y,t) = -1;
                end
                if(bb*HH<=0 || (b-pi)*(H-pi)<0)
                    lvl(x,y,t) = 1;
                    lvl(x,y-1,t) = 1;    
                elseif((b-hp)*(H-hp)<0 || (b-dp)*(H-dp)<0)
                    lvl(x,y,t) = -1;
                    lvl(x,y-1,t) = -1;
                end
            end
        end
    end
    close(w);
end

function lvl = makelvl2(Adiff, Adiff2, maskout)
    % for mode 2
    global nx ny nt
    lvl = zeros(nx,ny,nt,'single'); % level sets
    w = waitbar(0,'Building level sets...');
    for x=2:nx
        waitbar(x/(nx-1));
        for y=2:ny
            if maskout(x,y)
                continue
            end
            for t=1:nt
                a1 = Adiff(x-1,y,t);
                b1 = Adiff(x,y-1,t);
                H1 = Adiff(x,y,t);
                a2 = Adiff2(x-1,y,t);
                b2 = Adiff2(x,y-1,t);
                H2 = Adiff2(x,y,t);
                if(a1*H1<0)
                    lvl(x,y,t) = 1;
                    lvl(x-1,y,t) = 1;
                elseif(a2*H2<0)
                    lvl(x,y,t) = -1;
                    lvl(x-1,y,t) = -1;
                end
                if(b1*H1<0)
                    lvl(x,y,t) = 1;
                    lvl(x,y-1,t) = 1;    
                elseif(b2*H2<0)
                    lvl(x,y,t) = -1;
                    lvl(x,y-1,t) = -1;
                end
            end
        end
    end
    close(w);
end

function [D1, D2] = dist(lvl) 
    global nx ny nt
    D1 = zeros(nx,ny,nt,'single');
    D2 = D1;
    w = waitbar(0,'Computing distance functions...');
    for t=1:nt
        if mod(t,50)==0
           waitbar(t/nt);
        end
        bw1 = (lvl(:,:,t)==1);
        bw2 = (lvl(:,:,t)==-1);
        D1(:,:,t) = bwdist(bw1);
        D2(:,:,t) = bwdist(bw2);
    end
    D1(isinf(D1)) = NaN;
    D2(isinf(D2)) = NaN;
    close(w);
end

function mo2 = edgemask(maskout)
    mo2 = maskout;
    mo2(1,:) = 1;
    mo2(end,:) = 1;
    mo2(:, 1) = 1;
    mo2(:, end) = 1;
end

function [num, ch, xc1, yc1, xc2, yc2, yi, xi, chi, I, J] = findPS(D1new, D2new, maskdist, t, dis)
    global nx ny
    num = 0;
    ch = 0;
    [C1, ~] = contour(D1new(:,:,t),[0,0],'-w','Visible','off');
    [C2, ~] = contour(D2new(:,:,t),[0,0],'-k','Visible','off');
    
    % separate contours into their component pieces
    [xc1, yc1] = contoursep(C1);
    [xc2, yc2] = contoursep(C2);
    
    % gradients of distance functions
    [gd1x, gd1y] = gradient(D1new(:,:,t));
    [gd2x, gd2y] = gradient(D2new(:,:,t));
    % find phase singularities
    if (size(xc1,2)>1 && size(xc2,2)>1)
        % for some reason x and y are switched!?
        [xi,yi,I,J] = intersections(xc1,yc1,xc2,yc2);
        % throw out PSs that are NaNs
        keep = ~isnan(I) & ~isnan(J);
        xi = xi(keep);
        yi = yi(keep);
        I = I(keep);
        J = J(keep);
        % compute chirality for each singularity   
        ni = size(xi,1);
        chi = zeros(ni,1);
        for z=1:size(xi,1)
            % where is this intersection?
            xloc = xi(z);
            yloc = yi(z);
            % skip this ps if it is close to the edge of mask
            if xloc<=1 || yloc<=1 || xloc>=nx-1 || yloc>=ny-1 ...
                || maskdist(round(xloc),round(yloc))<dis
                continue
            end
            neux = floor(xloc):ceil(xloc); % need meshgrid locally only
            if size(neux,2)==1 % integer case
                neux = [xloc xloc+1];
            end
            neuy = floor(yloc):ceil(yloc);
            if size(neuy,2)==1 % integer case
                neuy = [yloc yloc+1];
            end
            [Xnr,Ynr] = meshgrid(neux,neuy); % for local interpolation
            % find gradients
            gdh1x = interp2(Xnr,Ynr,gd1x(neux,neuy),xloc,yloc);
            gdh1y = interp2(Xnr,Ynr,gd1y(neux,neuy),xloc,yloc); 
            gdh2x = interp2(Xnr,Ynr,gd2x(neux,neuy),xloc,yloc);
            gdh2y = interp2(Xnr,Ynr,gd2y(neux,neuy),xloc,yloc); 
            chi(z) = sign(gdh1x*gdh2y-gdh1y*gdh2x); % z-component of cross product
            if chi(z)>0 % CW rotation - positive chi - black
               num = num+1;
               ch = ch+1;
            elseif chi(z)<0 % CCW rotation - negative chi - white
               num = num+1;
               ch = ch-1;
            end
        end
        return
    else % no phase singularities for some reason
        num=0; ch=0; xi=[]; yi=[]; chi=[]; I=[]; J=[];
    end   
end