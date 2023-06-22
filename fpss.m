function [xs ys c N chrg] = fpss(pht, k, nx, ny)
    xs = [];
    ys = xs;
    c = xs;
    N = 0;
    chrg = 0;
    for x=1:nx-1
        for y=1:ny-1
            sqr = pht(x:x+1,y:y+1);
            [chi, dx, dy] = localize(sqr, k);
            if chi~=0
               xs = [xs x+dx];
               ys = [ys y+dy];
               c = [c chi];
               N = N+1;
               chrg = chrg+chi;
            end
        end
    end
end