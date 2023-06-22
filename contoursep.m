function [xc, yc] = contoursep(C)
    index = 1;
    xc = C(2,:); % x and y are switched when plotting
    yc = C(1,:);
    while size(C,2)>index
        len = C(2,index); % number of points in this part
        xc(index) = NaN; % new separator
        yc(index) = NaN;
        index = index+len+1;
    end
end