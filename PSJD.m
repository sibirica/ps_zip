function [] = PSJD(filename)
    global nx ny nt
    file = sprintf('%s.mat', filename);
    fileJD = sprintf('%sJD.mat', filename);
    load(file, 'A', 'nx', 'ny', 'nt', 'T');
    tau = round(0.1125*T);
    sigma = 2; % for Gaussian filter of D (2)
    R = 3; % ~ min separation
    ker = fspecial('gaussian', R,1);
    disp('Computing Jacobian determinant...');
    D = JD(A, tau);
    D(1,:,:)=0;
    D(end,:,:)=0;
    D(:,1,:)=0;
    D(:,end,:)=0;
    D = imgaussfilt(D, sigma);
    disp('Finding peaks...')
    for t=1:nt-tau
        % the way the threshold is used is weird
        Dt = squeeze(D(:,:,t));
        thres = max(abs(Dt(:)))/3;
        pks = FastPeakFind(Dt,thres*2^16/max(Dt(:)),ker,2,1);
        n = length(pks)/2;
        xs(1:n,t) = pks(1:2:end);
        ys(1:n,t) = pks(2:2:end);
        chi(1:n,t) = -1;
        pks2 = FastPeakFind(-Dt,thres*2^16/max(Dt(:)),ker,2,1);
        m = length(pks2)/2;
        xs(n+1:n+m,t) = pks2(1:2:end);
        ys(n+1:n+m,t) = pks2(2:2:end);
        chi(n+1:n+m,t) = 1;
    end
    sv = 10; % velocity scaling factor
    disp('Tracking PSs!')
    [tracks, ~, ~] = trackPS(xs,ys,chi,sv);
    trat = [];
    if numel(tracks)==0
       fail = 1;
    else
       fail = 0;
    end
    for t=1:nt
       if ~fail
           trat{t} = tracks(tracks(:,4)==t,:);
       end
    end
    save(fileJD, 'D', 'xs', 'ys', 'chi', 'tracks', 'trat');
end

function D = JD(V,tau)
    global nx ny nt
    dvdx = cat(1, V(2,:,:)-V(1,:,:), (V(3:end,:,:)-V(1:end-2,:,:))/2, V(end,:,:)-V(end-1,:,:));
    dvdy = cat(2, V(:,2,:)-V(:,1,:), (V(:,3:end,:)-V(:,1:end-2,:))/2, V(:,end,:)-V(:,end-1,:));
    D = dvdx(:,:,1:end-tau).*dvdy(:,:,1+tau:end)-dvdx(:,:,1+tau:end).*dvdy(:,:,1:end-tau);
end