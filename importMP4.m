%clear all;
name='Mechanism6';
extension='mp4';
scale=2;
system('rm -f /tmp/frame*');
system(sprintf('ffmpeg -i %s.%s -q:v 1 -vsync 0 -vf "scale=iw/%d:ih/%d" /tmp/frame%%05d.jpg',name,extension,scale,scale));
[status,cmdout] = system('ls -1 /tmp/frame* | wc -l');
N=str2num(cmdout);
u=im2single(imread(sprintf('/tmp/frame00001.jpg',i)));
[L1 L2]=size(rgb2gray(u));
data=single(zeros(L1,L2,N));
for i=1:N  
    u=im2single(imread(sprintf('/tmp/frame%05d.jpg',i)));
    data(:,:,i)=rgb2gray(u);
end
save(sprintf('%s.mat',name),'data');
