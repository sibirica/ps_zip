function xdiff = tdiff(x)
    xdiff = cat(3, x(:,:,2)-x(:,:,1), (x(:,:,3:end)-x(:,:,1:end-2))/2, x(:,:,end)-x(:,:,end-1));
end