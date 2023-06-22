function [] = plotByPieces(C, option)
    if isempty(C)
       return
    end
    nans = find(isnan(C(1,:)'));
    if isempty(nans)
        plot(C(2,1:end),C(1,1:end), option, 'LineWidth', 2); 
        %LineWidth 3 is better for figures
    else
        inter1 = cat(1, 0, nans); % contains 0 and then location of each NaN
        NaNinds = cat(1, inter1, size(C,2)); % 0, locs of NaNs, end
        ind = 1;
        while ind<size(NaNinds,1)
            st = NaNinds(ind)+1;
            en = NaNinds(ind+1);
            ind = ind+1;
            plot(C(2,st:en),C(1,st:en), option, 'LineWidth', 2);
        end
    end

end