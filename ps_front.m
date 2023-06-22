% curvature of wavefront = div(grad(D2new)/norm(grad(D2new)))

function Df = ps_front(D1new, D2new, maskout)
    global nx ny nt
    W = 4; % (half-)width in pixels of front/back
    %sig = [0.01 0.01 6];
    disp('Computing conduction velocity...')
    Df = diff(D2new, 1,3);
%     twind = 5;
%     Df = zeros(size(D2new));
%     w = waitbar(0,'Differentiating...'); 
%     for x=1:nx
%         waitbar(x/nx);
%         for y=1:ny
%             [Df(x,y,:), ~] = diff2(squeeze(D2new(x,y,:)), twind);
%         end
%     end
%     close(w);
    %Df = imgaussfilt3(Df, sig);
    Df = abs(Df);
    Df(:,:,nt) = Df(:,:,nt-1);
    Dg = zeros(size(D2new));
    for t=1:nt
        [grx, gry] = gradient(D2new(:,:,t));
        Dg(:,:,t) = sqrt(grx.^2+gry.^2);
    end
    Df = min(Df./Dg,5);
    for x=1:nx
        for y=1:ny
            if maskout(x,y)
               Df(x,y,:) = 0;
               continue 
            end
            st = 1;
            flag = 0;
            avg = 0;
            for t=2:nt
                if abs(D2new(x,y,t)) > W || D1new(x,y,t)<0
                   if flag
                      avg = mean(Df(x,y,st:t-1));
                      Df(x,y,st:t-1) = avg;
                      flag = 0;
                   end
                   Df(x,y,t) = avg;
                   st = t+1;
                else
                   flag = 1;
                end
            end
        end
    end
    %figure(1);
    %plot(squeeze(Df(64,64,:)));
    %hold on
    %plot(squeeze(Dg(64,64,:)));
end