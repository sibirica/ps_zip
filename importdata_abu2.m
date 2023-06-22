clear all;
m=4; % reduce data size by this factor
n=0.25; % starting index
dir='/home/roman/Projects/cardiac/flavio/abubu/pig_2d_2_spirals/';
prefix='';
format='%s%s%.2f.png';
[status,cmdout] = system(sprintf('ls -1 %s | wc -l',dir));
N=str2num(cmdout);
u=im2single(imread(sprintf(format,dir,prefix,n)));
uc=rgb2gray(u);
[L1 L2]=size(uc(1:m:end,1:m:end));
data=single(zeros(L1,L2,N));
h=waitbar(0,'Processing');
for i=1:N
    waitbar(i/N,h);
    u=im2single(imread(sprintf(format,dir,prefix,(i-1)*0.25+n)));
    uc=rgb2gray(u);
    data(:,:,i)=uc(1:m:end,1:m:end);
end
save('pig_2d_2spirals.mat','data');
