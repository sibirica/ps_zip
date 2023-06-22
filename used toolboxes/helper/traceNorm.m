function normalized = traceNorm(trace,varargin)
% normalized.m takes in a Vm/Ca trace and outputs a normalized pixel. If
% the trace is just a single pixel, just run this shit. If you're passing
% in the whole fucking matrix of pixels, then put the row and column number
% too as arguments.

if length(size(trace)) == 2;
    tmp = trace - min(trace);
    normalized = tmp ./ max(tmp);
else
    if isempty(varargin);
        error('give a fucking pixel')
    end
    row = varargin{1};
    col = varargin{2};
    tmp = squeeze(trace(row,col,:)) - min(squeeze(trace(row,col,:)));
    normalized = tmp ./ max(tmp);
end

