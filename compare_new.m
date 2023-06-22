function [tscore, dscore] = compare_new(file1, file2, cut)
    % file1 should be the reference recording
    % can cut beginning and end of recording, but seems unnecessary
    alpha = 0.3; % fraction of L that counts as close enough
    beta = 0.8; % fraction of T that counts as long enough for match 
    load(strcat(file1,'.mat'), 'nx', 'ny', 'nt', 'T');
    load(strcat(file1,'Y.mat'), 'tracks', 'trat', 'L');
    trat1 = trat;
    tracks1 = tracks;
    load(strcat(file2,'Y.mat'), 'tracks', 'trat');
    %load(strcat(file2,'.mat'), 'tracks', 'trat');
    trat2 = trat;
    tracks2 = tracks;
    df1 = 0;
    df2 = 0;
    if nargin>2
       df1 = size(find(tracks1(:,4)<1+cut | tracks1(:,4)>nt-cut),1);
       df2 = size(find(tracks2(:,4)<1+cut | tracks2(:,4)>nt-cut),1);
    else
       cut = 0;
    end
    md = alpha*L;
    mt = beta*T;
    matches = [];
    tmatch = 0;
    d = 0;
    ttotal = size(tracks1,1)+size(tracks2,1)-df1-df2;
    for t=1+cut:nt-cut
        matches = [matches; findMatch(trat1, trat2, md, t)];       
    end
    allPairs = unique(matches(:,1:2),'rows');
    for k=1:size(allPairs,1)
        pair = allPairs(k,:);
        pms = matches(matches(:,1) == pair(1) & matches(:,2) == pair(2),:);
        L = size(pms,1);
        if L<mt
            continue
        end
        tmatch = tmatch + L;
        d = d+sum(pms(:,3));
    end
    tscore = tmatch*2/ttotal;
    dscore = d/tmatch;
end

function [matches] = findMatch(trat1, trat2, md, t)
    frm1 = trat1{t};
    frm2 = trat2{t};
    L1 = size(frm1,1);
    L2 = size(frm2,1);
    matches = [];
    while L1>0
        ds = zeros(L2,1);
        x1 = frm1(1,1);
        y1 = frm1(1,2);
        for j=1:L2
            x2 = frm2(j,1);
            y2 = frm2(j,2);
            if frm1(1,3)-frm2(j,3)==0 % same chiralities
                ds(j)=sqrt((x1-x2)^2+(y1-y2)^2);
            else
                ds(j)=Inf; % shouldn't be matched
            end
        end
        [M,I] = min(ds);
        if M<md
            % last entry might turn out to be unnecessary?
            matches = [matches; frm1(1,5) frm2(I,5) M t]; %id1 id2 dist t
            frm1(1,:) = [];
            frm2(I,:) = [];
        else
            frm1(1,:) = [];
        end
        L1 = size(frm1,1);
        L2 = size(frm2,1);
    end
end