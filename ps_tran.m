function [cre, ane, typeCr, typeAn] = ps_tran(tracks, grps, vmax, nt)
    global T
    md = 2*vmax; %maximum distance between paired PSs
    [cre, ane] = cran(tracks,md,nt); % find all events
    lcr = size(cre,1);
    lan = size(ane,1);
    typeCr = zeros(lcr,1);
    typeAn = zeros(lan,1);
    for i=1:lcr
       typeCr(i)=classifyCr(cre(i,:),grps);
    end
    for i=1:lan
       typeAn(i)=classifyAn(ane(i,:),grps);
    end
end

function [cre, ane] = cran(tracks,md,nt)
    lm = 10; % minimum length of a track
    cre = [];
    ane = [];
    cr = cell(nt,1);
    an = cell(nt,1);
    ma = max(tracks(:,5)); % number of ids
    for id=1:ma
        st = tracks(find(tracks(:,5)==id,1,'first'),:);
        en = tracks(find(tracks(:,5)==id,1,'last'),:);
        t1 = st(4);
        t2 = en(4);
        if t2-t1>=lm
            cr{t1} = [cr{t1}; st];
            an{t2} = [an{t2}; en];
        end
    end
    for t=2:nt-1
        ec = pair(cr{t},md);
        for i=1:size(ec)
            cre = [cre; ec(i,:) t];
        end
        ea = pair(an{t},md);
        for i=1:size(ea)
            ane = [ane; ea(i,:) t]; 
        end
    end
end

function [e]=pair(ps,md)
    % each row of e contains 2 ids of paired PSs
    lc = size(ps,1);
    e = [];
    while lc>1
        ds = zeros(lc-1,1);
        for j=2:lc
            x1 = ps(1,1);
            y1 = ps(1,2);
            x2 = ps(j,1);
            y2 = ps(j,2);
            if ps(1,3)+ps(j,3)==0 % opposite chiralities
                ds(j-1)=sqrt((x1-x2)^2+(y1-y2)^2);
            else
                ds(j-1)=Inf; % shouldn't be matched
            end
        end
        [M,I] = min(ds);
        if M<md
           e = [e; ps(1,5) ps(I+1,5)];
           ps(1,:) = [];
           ps(I,:) = [];
        else
           ps(1,:) = [];
        end
        lc = size(ps,1);
    end
end

function [type]=classifyAn(e,grps)
    ell = 5;
    x = [0 0 0 0];
    while min(x(:,4))<1
        frm = grps{e(3)};
        if isempty(frm)
            break
        end
        e(3) = e(3)-1;
        x = frm((frm(:,2)==e(1) & frm(:,3)==e(2)) | (frm(:,2)==e(2) & frm(:,3)==e(1)),:);
        x(x(:,5)./x(:,4)>ell, :)=[];
        if isempty(x)
            break
        end
    end
    a = any(x(:,1)==1);
    b = any(x(:,1)==2);
    c = any(x(:,1)==3);
    d = any(x(:,1)==4);
    if a+b+c+d<2
       fprintf("Weird annihilation: 0 or 1 level sets at t=%d \n",e(3)+1);
       fprintf("a: %d, b: %d, c: %d, d: %d \n",a,b,c,d);
       type = 0;
       return
    end
    if ~a
       if ~c
           type = 2; % b
       elseif ~d
           type = 8; % h
       else
           type = 1; % a
       end
    elseif ~b
       if ~c
           type = 4; % d
       elseif ~d
           type = 6; % f
       else
           type = 5; % e
       end
    elseif ~c
        type = 3; % c
    elseif ~d
        type = 7; % g
    else
        type = 9; % doubly degenerate
    end
    % f, g, h - merger
    % otherwise - collapse
end

function [type]=classifyCr(e,grps)
    global T
    ell = 5; % ratio of ad to dist that is "long"
    x = [0 0 0 0];
    while min(x(:,4))<1
        frm = grps{e(3)};
        if isempty(frm)
            break
        end
        e(3) = e(3)+1;
        x = frm((frm(:,2)==e(1) & frm(:,3)==e(2)) | (frm(:,2)==e(2) & frm(:,3)==e(1)),:);
        x(x(:,5)./x(:,4)>ell, :)=[];
        if isempty(x)
            break
        end
    end
    a = any(x(:,1)==1);
    b = any(x(:,1)==2);
    c = any(x(:,1)==3);
    d = any(x(:,1)==4);
    if a+b+c+d>2
       fprintf("Weird creation: 3 or 4 level sets at t=%d \n",e(3)-1);
       fprintf("a: %d, b: %d, c: %d, d: %d \n",a,b,c,d);
       type = 0;
       return
    end
    if a
       if c
           type = 6; % f
       elseif d
           type = 4; % d
       else
           type = 5; % e
       end
    elseif b
       if c
           type = 8; % h
       elseif d
           type = 2; % b
       else
           type = 1; % a
       end
    elseif c
        type = 7; % g
    elseif d
        type = 3; % c
    else
        type = 9; % doubly degen
    end
%     if a && ~d
%         type = 1;
%     elseif b && ~c
%         type = 2;
%     elseif c
%         type = 3;
%     else
%         type = 4;
%     end
    % 1 - wave coalescence
    % 2 - wave breakup 
    % 3 - conduction block
    % 4 - back ignition
end