function [ movies, headers, times ] = opread_g( file, cameraIndices, machineformat)
%OPREAD_G reads a movie file, version 'g'
%   Arguments are:
%   file            The filename
%   cameraIndices   The indices of the uses cameras, defaults to
%                   [1 2 3 ...].
%   machineformat   The machineformat (byte order), defaults to 'l' (little
%                   endian).
% 
%   Returns a cell array with the movies and optionally the headers and a
%   cell array with the exposure times.

if (nargin == 0)
  error('Too few arguments.');
end

if (nargin < 3)
  machineformat = 'l'; % default format = little endian
  if (nargin < 2)
    cameraIndices = 0:255;
  end
elseif (nargin > 3)
  disp(nargin);
  error('Too many arguments.');
end

% read headers
[cameraNumbers, headers] ...
  = opheader_g(file, cameraIndices, machineformat);
nCams = numel(cameraNumbers);

% read data

% open file (again)
[fid,msg] = fopen(file, 'rb', machineformat);
if fid == -1
   error('Cannot open file %s -> %s', file, msg);
end

fseek(fid, 9 + (1024+4)*nCams, 'bof'); % skip headers

% preallocate memory:
times = cell(1, nCams);
movies = cell(1, nCams);
for i = 1:nCams
  times{i} = zeros(1, headers{i}.frames, 'uint64');
  movies{i} = zeros(headers{i}.sizeY, ...
    headers{i}.sizeX, ...
    headers{i}.frames, ...
    'uint16');
end

%loop over whole data file
frames = zeros(1, nCams);
while true
  % from which camera is the next image?
  camIdx = fread(fid, 1, '*uint32');
  if (feof(fid))
    disp('EOF reached.');
    break;
  end
  internalIdx = find(cameraNumbers == camIdx, 1);
  
  frames(internalIdx) = frames(internalIdx) + 1;
  h = headers{internalIdx};
  
	image = fread(fid, [h.sizeX, h.sizeY], '*uint16'); 
	image = image'; % Data comes out transposed.
	image(1,1) = 65535;
  % Write the result as the next frame in the data matrix.
  movies{internalIdx}(:,:, frames(internalIdx)) = cast(image, 'uint16');
	if h.version == 'f'
		times{internalIdx}(frames(internalIdx)) = fread(fid, 1, '*uint64');
	end
% 	if (mod(loop,50) == 0)
% 		stop = progressbar(loop/(options.toFrame-options.fromFrame+1));
% 		if (stop) 
% 			disp('OPREAD aborted by user.');
% 			break;
% 		end
% 	end
end
% progressbar(1); % finalize the progress bar.
% close(gcf); % close the progressbar figure.
fclose(fid);

end

