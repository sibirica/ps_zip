function out = smoothers(data, varargin)

switch varargin{1}
    case 'time';
        disp('Time smoothing...');
        
        %% time smooth
        if length(varargin) == 1;
            windowSize = 1;
        else
            windowSize = varargin{2};
        end
        a = 1;
        b = (1./windowSize).*ones(1,windowSize);

        out = data;

        for row = 1:128
            for col = 1:128
                out(row, col, :) = ...
                    filter(b, a, squeeze(data(row, col, :)));
            end
        end
        out(:, :, 1) = out(:, :, 2);

    case 'space';
        disp('Space smoothing...');

        %% space smooth

        if length(varargin) >= 3;
            r = varargin{2};
            sigma = varargin{3};
        else
            r = 3;
            sigma = 3;
        end
        spatialFilter = fspecial('gaussian', [r r], sigma);

        out = imfilter(data, spatialFilter, 'replicate');

    otherwise
        disp('not acceptable smoother')
end

end