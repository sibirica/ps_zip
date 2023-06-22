function [] = findPSsparse(filename)
    file = sprintf("%s.mat",filename);
    fileX = sprintf("%sX.mat",filename);
    disp("Loading saved data...")
    load(file, 'ph', 'k', 'nx', 'ny', 'nt', 'sc');
    
    numPS = zeros(nt,1);
    charge = numPS;
    xs = [];
    ys = xs;
    chi = xs;
    w = waitbar(0,'Finding phase singularities...');
    for t=1:nt
        if mod(t,50)==0
           waitbar(t/nt); 
        end
        pht = ph(:,:,t);
        % analysis for only one frame on scaled grid
        [x y c N chrg] = fpss(pht, k*2^sc, nx/2^sc+1, ny/2^sc+1);
        xs(1:length(x),t) = (x-1)*2^sc+1;
        ys(1:length(x),t) = (y-1)*2^sc+1;
        chi(1:length(x),t) = c;
        numPS(t) = N;
        charge(t) = chrg;
    end
    close(w);
    
    disp("Saving PSs...")
    save(fileX, 'numPS', 'charge', 'xs', 'ys', 'chi');
end