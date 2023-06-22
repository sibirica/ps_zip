function [As A] = sparsify(A, sc)
    global nx ny nt
    A(nx+1,:,:) = A(nx,:,:);
    A(:,ny+1,:) = A(:,ny,:);
    As = A;
    A = A(1:2^sc:end,1:2^sc:end,:);
    w = waitbar(0,'Sparsifying...');
    for t=1:nt
        if mod(t,50)==0
           waitbar(t/nt);
        end
        %As(:,:,t) = interp2(A(:,:,t), sc, 'makima'); % used spline for paper, I think
        As(:,:,t) = interp2(A(:,:,t), sc, 'bicubic');
    end
    close(w);
    As = As(1:nx,1:ny,:);
end