function [C1, C2, grp] = repackage(xc, yc, D, tr, xp, yp, inds, sh, mode)
% 1) separates the points of the contour given by (xc,yc) into the two
% level sets based on distance function.
% 2) computes the graph associated with each level set
% row of grp: [type id1 id2 dist arcdist]
% Note: sh is shift in type based on I or J
    en = find((xc(:)==0) & (yc(:)==0),1);
    if ~isempty(en)
        xc = xc(1:en-1); % get rid of the the trailing 0s
        yc = yc(1:en-1);
    end
    si = size(xc,1);
    if si==0
        C1 = [];
        C2 = [];
        grp = [];
        return
    end
    C1 = zeros(2,ceil(si*2/3));
    C2 = C1;
    naninds = []; % where are the NaNs
    ind1 = 1; % index at which we are inserting into the contour
    ind2 = 1;
    ind3 = 1; % same for naninds
    prev = 0; % last level set added to
    for i=1:si
        if isnan(xc(i))
           C1(:,ind1)=[NaN NaN];
           C2(:,ind2)=[NaN NaN];
           naninds(ind3) = i;
           ind1 = ind1+1;
           ind2 = ind2+1;
           ind3 = ind3+1;
           prev = 0;
        elseif xc(i)<=0 || yc(i)<=0
            continue
        % D2, cos(ph) is negative for pi, D1, sin(ph) for 3pi/2
        elseif D(round(xc(i)), round(yc(i)))<0 % (x and y not switched)
           if prev==1
              C1(:,ind1)=[NaN NaN];
              ind1 = ind1+1;
           end
           C2(:,ind2) = [xc(i) yc(i)];
           ind2 = ind2+1;
           prev = 2;
        else
           if prev==2
              C2(:,ind2)=[NaN NaN];
              ind2 = ind2+1;
           end
           C1(:,ind1) = [xc(i) yc(i)];
           ind1 = ind1+1;
           prev = 1;
        end
    end
    % cropping out zeros
    C1 = C1(:,1:ind1-1);
    C2 = C2(:,1:ind2-1);
    
    % 2)
    dmax = 2; % how much skipped distance allowed
    grp = [];
    if numel(tr)==0 % tracking failed
        return
    end
    en = find(inds==0,1);
    if ~isempty(en)
        inds = inds(1:en-1);
    end
    naninds = [naninds(2:end)]';
    s = size(naninds,1);
    [Q,I] = sort(inds);
    lq = size(Q,1);
    x = [];
    y = [];
    TRl = [];
    for j=1:lq
    	Bj = I(j);
        x(j) = xp(Bj);
        y(j) = yp(Bj);
        trj = find(tr(:,1)==x(j) & tr(:,2)==y(j));
        if isempty(trj)
           trj=0;
        end
        TRl(j) = trj;
    end
    % remove PSs that dropped out
    kp = TRl(:)~=0;
    Q = Q(kp);
    x = x(kp);
    y = y(kp);
    TRl = TRl(kp);
    lq = size(Q,1);
    if lq < 2
       return 
    end
    iN = 1; % where we are in the NaNs
    N = 0;
    iG = 1; % where we are in grp
    TRs = 0;
    TRe = 0;
    for indQ=1:lq-1
        flag = 1; % 0 if there is definitely no edge
        Q1=Q(indQ);
        Q2=Q(indQ+1);
        x1=x(indQ);
        y1=y(indQ);
        x2=x(indQ+1);
        y2=y(indQ+1);
        TR1 = TRl(indQ);
        TR2 = TRl(indQ+1);
        if s>0
            N = naninds(iN);
        end
        ex = 0; % whether we are ready to exit
        while N<Q2 && ex==0 && s>0
            if N>Q1
               if TRs~=0 && flag==1 % normal loop
                  d = sqrt((xc(N-1)-xs)^2+(yc(N-1)-ys)^2); % leap when we loop
                  if d<dmax && TR1~=TRs
                    id1 = tr(TRs,5);
                    id2 = tr(TR1,5);
                    dist = sqrt((tr(TRs,1)-x1)^2+(tr(TRs,2)-y1)^2);
                    ad = (Qs-Qn)+(N-Q1)+d; % piece 1 + piece 2 + leap handicap
                    md = (D(round(xs),round(ys))<0);
                    type = 1+sh*2+md; % corresponds to C1-C4; (+4)
                    grp(iG,:) = [type id1 id2 dist ad];
                    iG = iG+1; 
                  end
               elseif iN==1
                   TRe=TR1;
                   Qe=Q1;
                   Qm=N-1;
                   xe=xc(i-1);
                   ye=yc(i-1);
               end 
               TRs=TR2; % reference to starting TR
               Qs=Q2;
               Qn=N+1; % index of the old NaN
               xs=xc(Qn); % position next to NaN, not at TRs
               ys=yc(Qn);
               flag = 0;
               iN = min(iN+1,s);
               N = naninds(iN);
            else
               iN = min(iN+1,s);
               N = naninds(iN);
            end
            if iN==s
                ex=1;
            end
        end
        if flag==1
           id1 = tr(TR1,5);
           id2 = tr(TR2,5);
           dist = sqrt((x1-x2)^2+(y1-y2)^2);
           ad = Q2-Q1;
           rp = round((Q1+Q2)/2);
           if isnan(xc(rp))
               rp=rp+1;
           end
           df = D(round(xc(rp)),round(yc(rp)));
           if abs(df)<1 % too close to tell
               a = 0;
               b = 0;
               p1 = floor(Q1)-1;
               p2 = ceil(Q2)+1;
               if p1>=2 && ~isnan(xc(p1))
                   a = -D(round(xc(p1)),round(yc(p1)));
               end
               if p2<=si && ~isnan(xc(p2))
                   b = -D(round(xc(p2)),round(yc(p2)));
               end
               df = (a+b)/2;
           end
           md = (df<0);
           type = 1+sh*2+md; % corresponds to C1-C4;
           grp(iG,:) = [type id1 id2 dist ad];
           iG = iG+1;
        end
    end
    % edge case: looping around from the end to the beginning
    if s>0
        Na = naninds(iN);
    else
        Na = 0;
    end
    la = 0; % looparound flag
    if TRe~=0 
       d = sqrt((xc(end)-xc(2))^2+(yc(end)-yc(2))^2); % leap when we loop
       if Na<Q2 && d<dmax
           id1 = tr(TRs,5);
           id2 = tr(TRe,5);
           dist = sqrt((tr(TRs,1)-tr(TRe,1))^2+(tr(TRs,2)-tr(TRe,2))^2);
           ad = si-(Qs-Qe)+d;
           md = (D(round(xe),round(ye))<0);
           type = 1+sh*2+md; % corresponds to C1-C4; (+4)
           grp(iG,:) = [type id1 id2 dist ad];
           la = 1;
       elseif TRe~=TRl(1)
           d = sqrt((xe-xc(2))^2+(ye-yc(2))^2);
           if d<dmax
               id1 = tr(TRl(1),5);
               id2 = tr(TRe,5);
               dist = sqrt((tr(TRl(1),1)-tr(TRe,1))^2+(tr(TRl(1),2)-tr(TRe,2))^2);
               ad = (Q(1)-2)+(Qm-Qe)+d;
               md = (D(round(xe),round(ye))<0);
               type = 1+sh*2+md; % corresponds to C1-C4; (+4)
               grp(iG,:) = [type id1 id2 dist ad];
               iG = iG+1;
           end
       end
    end
    % edge case: loop includes last PS
    if Na<Q2 % iN==s || s==0 
        N = si-1; % position before ending NaN
    else
        %N = naninds(iN+1);
        N = naninds(find(naninds(:)>Q2,1))-1; % position before ending NaN
    end
    xe = xc(N);
    ye = yc(N);
    if TRs==0
       TRs=TRl(1);
       Qs=Q(1);
    end
    iM = find(naninds(:)<Qs,1,'last');
    if isempty(iM)
        M = 2;
    else
        M = naninds(iM)+1; % position after starting NaN
    end
    xs = xc(M);
    ys = yc(M);
    d = sqrt((xe-xs)^2+(ye-ys)^2); % leap when we loop
    if d<dmax && TRs~=TR2 && ~la
        id1 = tr(TRs,5);
        id2 = tr(TR2,5);
        dist = sqrt((tr(TRs,1)-x2)^2+(tr(TRs,2)-y2)^2);
        ad = (Qs-M)+(N-Q2)+d;
        md = (D(round(xs),round(ys))<0);
        type = 1+sh*2+md; % corresponds to C1-C4; (+4)
        grp(iG,:) = [type id1 id2 dist ad];
    end
    % fix duplicates (dup glitch) by flipping type
    [~,ia,~]=unique(grp,'rows','stable');
    if ~isempty(ia)
        grp(ia,1)=bitxor(grp(ia,1)+1,1)-1;
    end
    if ~isempty(grp) && mode==2 % backwards for some reason
        grp(:,1)=bitxor(grp(:,1)+1,1)-1;
    end
end