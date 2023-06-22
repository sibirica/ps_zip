function [pathx, pathy, patht, len] = path2(tracks, id)
% follow trajectory of particle with certain id
inds = find(tracks(:,5)==id);
be = min(inds);
en = max(inds);
pathx = tracks(be:en,1);
pathy = tracks(be:en,2);
patht = tracks(be:en,4);
len = en-be+1;