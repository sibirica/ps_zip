i0 = 3800;
T = 4150;
L = T-i0;
m=4; % reduce data size by this factor
dir='PS/Claire/LEAP simulation/TenTusscher_2D_LEAPfib_L20/data/';
prefix='00';
format='%s%s%d_0(1).png';
%str = sprintf(format,dir,prefix,i0)
u=im2single(imread(sprintf(format,dir,prefix,i0)));
%uc=rgb2gray(u);
[L1, L2]=size(u(1:m:end,1:m:end));
data=single(zeros(L1,L2,L));
h=waitbar(0,'Processing');
for i=i0:T
    waitbar(i/L,h);
    u=im2single(imread(sprintf(format,dir,prefix,i)));
    %uc=rgb2gray(u);
    data(:,:,i)=u(1:m:end,1:m:end);
end
save('data.mat', 'data');