function [dy, dy2] = diff2(y, N)
% y - to be differentiated
% N - number of points used
    s = size(y);
    S = max(s);
    if s(1)<s(2)
        y = y';
        s = size(y);
    end
    dy = zeros(s);
    dy2 = dy;
    h = floor(N/2);
    r = (1:N)';
    for i=1:h
        [dy(i), dy2(i)] = fitDiff(r, y(r), i-h-1, 1);
    end
    for i=h+1:S-h
        r = (i-h:i+h)';
        [dy(i), dy2(i)] = fitDiff(r, y(r), 0, 1); 
    end
    r = (S-N+1:S)';
    for i=S-h+1:S
        [dy(i), dy2(i)] = fitDiff(r, y(r), i-S+h, 1); 
    end
end