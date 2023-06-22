% Return the data in a .dat file as a three-dimensional array 
% (image sequence), where the third coordinate is the frame number.
%
% v = OPREAD(file) returns the whole .dat file as a uint16 matrix.
% v = OPREAD(file, options) uses the structure options
%       .fromFrame     - number of the first frame to read
%       .toFrame       - number of the last frame to read
%       .fromX         - number of the first column to read
%       .toX           - number of the last column to read
%       .fromY         - number of the first row to read
%       .toY           - number of the last row to read
%       .machineformat - machineformat to be used for reading from the file
%       .precision     - returns the matrix as array of class 'precision'
% 
% Pixel and frame numbers can be specified with normal positive indices
% starting at 1, or with negative indices starting at 0 which corresponds
% to the greatest posssible index.
%
% [v,h] = OPREAD(...) returns also the header structure obtained from
%                     OPHEADER().
%
% [v,h,times] = OPREAD(...) returns also the times of exposure for each 
%                           frame as a cell-array of strings (versions E)
%                           or as a row vector of type uint64 (version F).
%
% Based on files written by Martin, Cornell University, 2005.
% Various modifications by Stefan Luther and Amgad Squires, 2006-2009.
% Unification and version F by Johannes Schroeder-Schetelig, 2010-01-12.

function [v,h,times] = opread(file, options)

if (nargin<1)
    error('Too few arguments. OPREAD requires at least a file name.');
end
if nargin > 2
	error('Too much arguments.');
end

if (nargin < 2)
	options = struct; % 1-by-1 structure with no fields.
end

if (~isfield(options, 'fromFrame'))
	options.fromFrame = 1;
end
if (~isfield(options, 'toFrame'))
	options.toFrame = 0;
end
if (~isfield(options, 'fromX'))
	options.fromX = 1;
end
if (~isfield(options, 'toX'))
	options.toX = 0;
end
if (~isfield(options, 'fromY'))
	options.fromY = 1;
end
if (~isfield(options, 'toY'))
	options.toY = 0;
end
if (~isfield(options, 'machineformat'))
	options.machineformat = 'l';
end
if (~isfield(options, 'precision'))
	options.precision = 'uint16';
end


h = opheader(file, options.machineformat);
disp('File header:');
disp(h);

if (h.version~='a' && h.version~='c' && h.version~='d' && h.version~='e' && h.version~='f')
    error('"%c" is not a recognized version.',h.version);
end


if (options.fromFrame <= 0)
	options.fromFrame = h.frames + options.fromFrame;
end
if (options.toFrame <= 0)
	options.toFrame = h.frames + options.toFrame;
end
if (options.fromX <= 0)
	options.fromX = h.sizeX + options.fromX;
end
if (options.toX <= 0)
	options.toX = h.sizeX + options.toX;
end
if (options.fromY <= 0)
	options.fromY = h.sizeY + options.fromY;
end
if (options.toY <= 0)
	options.toY = h.sizeY + options.toY;
end

% Bounding checking:
if (options.fromFrame < 1 || options.fromFrame > h.frames)
	error('options.fromFrame is out of bounds [1, %d].',h.frames);
end
if (options.toFrame < 1 || options.toFrame > h.frames)
	error('options.toFrame is out of bounds [1, %d].',h.frames);
end
if (options.fromX < 1 || options.fromX > h.sizeX)
	error('options.fromX is out of bounds [1, %d].',h.sizeX);
end
if (options.toX < 1 || options.toX > h.sizeX)
	error('options.toX is out of bounds [1, %d].',h.sizeX);
end
if (options.fromY < 1 || options.fromY > h.sizeY)
	error('options.fromY is out of bounds [1, %d].',h.sizeY);
end
if (options.toY < 1 || options.toY > h.sizeY)
	error('options.toY is out of bounds [1, %d].',h.sizeY);
end

disp(sprintf('Reading video range [Y=%d:%d, X=%d:%d, frame=%d:%d]', ...
	options.fromY, options.toY, options.fromX, options.toX, ...
	options.fromFrame, options.toFrame));


[fid,msg] = fopen(file,'rb',h.machineformat);
if fid == -1
   error('Cannot open file %s -> %s',file,msg);
end

fseek(fid, 1024, 'bof'); % skip header

if (h.version == 'e')
	fseek(fid, 6*options.fromFrame-1, 'cof');
	times = {};
	for i = 1:options.toFrame-options.fromFrame+1
		str = transpose(fread(fid, 6, 'int8=>char'));
	    times = [times str];
	end
	fseek(fid, 1024+6*h.frames, 'bof');
end

% preallocate memory:

if (h.version == 'f')
	times = zeros(1, options.toFrame-options.fromFrame+1, 'uint64');
end

v = zeros(options.toY-options.fromY+1, options.toX-options.fromX+1, ...
		  options.toFrame-options.fromFrame+1, options.precision);

		  
	  
fseek(fid, h.blockSize*(options.fromFrame-1), 'cof'); % skip first blocks
for loop = 1:(options.toFrame-options.fromFrame+1)         
	image = fread(fid, [h.sizeX, h.sizeY], '*uint16'); 
	image = image'; % Data comes out transposed.
	image(1,1) = 65535;
	image = image(options.fromY:options.toY, options.fromX:options.toX); % Trim the image to the desired dimensions.
	v(:,:,loop) = cast(image, options.precision); % Write the result as the next frame in the data matrix.
	if h.version == 'f'
		times(loop) = fread(fid, 1, '*uint64');
	end
	if (mod(loop,50) == 0)
		stop = 0;
		if (stop) 
			disp('OPREAD aborted by user.');
			break; 
		end
	end
end
fclose(fid);
return;
