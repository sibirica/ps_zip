function [] = run_batch()

load('karma.mat');

%opts.exp = 0;
%opts.noise = 0;
%opts.sc = 0;
%disp('--- Benchmark ---')
%ps_main(opts,'KarmaB',data,'mask');

opts.exp = 2;
opts.record = 0;
ns = [0 0.1 0.3 1];
sc = [0 2 3 4 5];
for i=1:size(sc,2)
    opts.sc = sc(i);
    for j=1:size(ns,2)
        %filename = strcat('KarmaN',num2str(j),'S',num2str(i));
        filename = strcat('LarmaN',num2str(j),'S',num2str(i));
        file = strcat(filename,'.mat');
        if exist(file,'file')
           continue; 
        end
        opts.noise = ns(j);
        resStr = num2str(size(data,1)/2^sc(i));
        disp(strcat(['--- Noise: '  num2str(ns(j))  ' ---']));
        disp(strcat(['--- Resolution: ' resStr 'x' resStr ' ---']));
        %fail = ps_main(opts,strcat('KarmaN',num2str(j),'S',num2str(i)),data,'mask');
        fail = ps_main(opts,strcat('LarmaN',num2str(j),'S',num2str(i)),data,'mask');
        if fail % noise too high to continue
           break; 
        end
    end
end