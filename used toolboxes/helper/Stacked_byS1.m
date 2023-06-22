function Stacked = Stacked_byS1(DATA,S1,pacing_interval)
% S1: the timing of the each pacing beat. It is a vector containing the
% time of each upstroke.
%
% pacing_interval: is the interval between two pacing beats.

if ndims(DATA) == 2;
    Stacked = zeros(1, pacing_interval);
    for i = 1:length(S1)
    Stacked = Stacked + DATA( S1(i) : S1(i) + pacing_interval - 1 )';%-mean(DATA( S1(i):S1(i)+pacing-1)');
    end
end

if ndims(DATA) == 3
    Stacked = zeros(size(DATA, 1), size(DATA, 2), pacing_interval);
    for i = 1:length(S1)
        Stacked(:, :, 1:pacing_interval) = Stacked(:, :, 1:pacing_interval)...
            + DATA(:, :, S1(i) : S1(i) + pacing_interval - 1);
        %-repmat( mean(DATA(:,:,S1(i):S1(i)+pacing-1),3), [1 ,1, pacing]);
    end
end

Stacked = Stacked / (length(S1) - 1);