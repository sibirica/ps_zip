function s = smush(xs, ys, chi)
% converts xs, ys, chi format to smushed 2d matrix for track
    [mPS, nt] = size(xs); % max # PS, number of time frames
    ind = 1;
    s = zeros(4,mPS/2*nt);
    for t=1:nt
        for n=1:mPS
            if xs(n,t)==0
                break
            elseif chi(n,t)~=0
                s(:,ind) = [xs(n,t),ys(n,t),chi(n,t)*100,t];
                ind = ind+1;
            end
        end
    end
    s = s(:,1:ind-1)';
end