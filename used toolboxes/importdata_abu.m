clear all;
m=4; % reduce data size by this factor
n=12000; % starting index
dir='/home/roman/Projects/cardiac/flavio/abubu/set_05/';
name='set_05_';
format='%s%s%5d.png';
[status,cmdout] = system(sprintf('ls -1 %s | wc -l',dir));
N=str2num(cmdout);
u=im2single(imread(sprintf(format,dir,name,n)));
uc=rgb2gray(u);
[L1 L2]=size(uc(1:m:end,1:m:end));
data=single(zeros(L1,L2,N));
for i=1:N  
    u=im2single(imread(sprintf(format,dir,name,i-1+n)));
    uc=rgb2gray(u);
    data(:,:,i)=uc(1:m:end,1:m:end);
end
save(sprintf('%s.mat',name),'data');
