function [chi, dx, dy] = localize(sqr, k)
% square is [0,1]x[0,1]
    %options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',0.01);
    options = optimset('TolX',0.01);
    f = zeros(1,2);
    v = zeros(3,2);
    for i=1:2
        fun = @(vars) err(vars(1), vars(2), vars(3), sqr, k, 3-2*i);
        xl = -2:0.4:3;
        yl = -2:0.4:3;
        phil = 0:1:2*pi;
        La = length(xl);
        Lb = length(yl);
        Lc = length(phil);
        errs = zeros(La,Lb,Lc);
        for a=1:La
            for b=1:Lb
                for c=1:Lc
                    errs(a,b,c)=fun([xl(a) yl(b) phil(c)]);
                end
            end
        end
        [M I] = min(errs(:));
        [s1 s2 s3] = ind2sub([La Lb Lc],I);
        x0 = [xl(s1) yl(s2) phil(s3)];
        [v(:,i) f(i)] = fminsearch(fun,x0,options); %[-0.5 0.5 pi]
    end
    if f(1)<f(2)
        chi = 1;
        dx = v(1,1);
        dy = v(2,1);
    else
        chi = -1;
        dx = v(1,2);
        dy = v(2,2);
    end
    if abs(dx-0.5)>0.5 || abs(dy-0.5)>0.5
       chi = 0; 
    end
end

function [e] = err(dx,dy,phi0,sqr,k,chi)
    xj = [0 1; 0 1];
    yj = [0 0; 1 1];
    rj = (xj-dx).^2+(yj-dy).^2;
    thj = atan2(yj-dy, xj-dx);
    phi = phi0-k*rj+chi*thj;
    d = ds(sqr,phi);
    e = double(sum(d(:)));
end

function [d] = ds(phi1, phi2)
    %d = cat(3,(phi1-phi2).^2,(phi1-phi2+2*pi).^2,(phi1-phi2-2*pi).^2);
    %d = cat(3,abs(phi1-phi2),abs(phi1-phi2+2*pi),abs(phi1-phi2-2*pi));
    %d = min(d,[],3);
    d = (phi1-mod(phi2,2*pi)).^2;
end