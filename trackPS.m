function [tr, L, vmax] = trackPS(xs, ys, chi, sv, Tin)
    % sv - velocity scaling factor
    global nx ny T
    if nargin<5
       Tin = T;
    end
    L = sqrt(nx*ny/size(xs,1));
    s = smush(xs,ys,chi);
    param.mem = 1;
    param.dim = 2;
    param.good = floor(Tin/4); % T/4 seems good for pig
    param.quiet = 1;
    vmax = sv*L/Tin;
    tr = [];
    if size(xs,1)>1 % otherwise, nothing to track
        tr = track(s,vmax,param); % format: x, y, 100*chi, t, id
    end
end