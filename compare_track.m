function [score, match, total, actual] = compare_track(file1, file2, alpha)
    % file1 should generally be the reference recording
    % L ~ the avg PS separation
    % alpha - fraction of L that counts as close enough
    % alpha = 0 - calculates avg distance between PSs in file1 and file2
    load(strcat(file1,'.mat'), 'nx', 'ny');
    load(strcat(file1,'Y.mat'), 'trat', 'L');
    trat1 = trat;
    load(strcat(file2,'Y.mat'), 'trat');
    trat2 = trat;
    if alpha~=0
        md = alpha*L;
        nt = size(trat1,1);
        match = zeros(nt,1);
        total = match;
        for t=1:nt
            frm1 = trat1{t};
            frm2 = trat2{t};
            L1 = size(frm1,1);
            actual(t) = L1;
            L2 = size(frm2,1);
            total(t) = L1+L2;
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
                   frm1(1,:) = [];
                   frm2(I,:) = [];
                   match(t) = match(t)+2;
                else
                   frm1(1,:) = [];
                end
                L1 = size(frm1,1);
                L2 = size(frm2,1);
            end
        end
        score = sum(match)/sum(total);
    else
        nt = size(trat1,1);
        match = zeros(nt,1);
        total = match;
        actual = total;
        for t=1:nt
            frm1 = trat1{t};
            frm2 = trat2{t};
            L1 = size(frm1,1);
            actual(t) = L1;
            L2 = size(frm2,1);
            if L2>L1
                big = frm2;
                small = frm1;
                L1 = size(frm2,1);
                L2 = size(frm1,1);
            else
                big = frm1;
                small = frm2;
            end
            total(t) = L1;
            if L2==0
               continue 
            end
            for i=1:L1
                ds = zeros(L2,1);
                x1 = big(i,1);
                y1 = big(i,2);
                for j=1:L2
                    x2 = small(j,1);
                    y2 = small(j,2);
                    if big(i,3)-small(j,3)==0 % same chiralities
                        ds(j)=sqrt((x1-x2)^2+(y1-y2)^2);
                    else
                        ds(j)=min([x1,y1,nx-x1,ny-y1]); % shouldn't be matched
                    end
                end
                M = min(ds);
                match(t) = match(t)+M;
            end
        end
        score = sum(match)/sum(total);
    end
end