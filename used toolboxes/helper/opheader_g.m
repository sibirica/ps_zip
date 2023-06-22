function [cameraNumbers, headers] = opheader_g(file, cameraIndices, machineformat)

% argument checking
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

% open the file
[fid,msg] = fopen(file,'rb',machineformat);
if fid<0
   error('Cannot open file %s -> %s',file,msg);
end

% read the headers
[cameraNumbers, headers] = read_g_header(fid, cameraIndices);

fclose(fid);

end

% returns an array of structs which have the same content as version 'f'
function [cameraNumbers, headers] = read_g_header(fid, cameraIndices)
version = fread(fid,1, 'int8=>char');
if (version ~= 'g')
  error('not version "g"')
end

byteOrderMark = fread(fid, 1, '*uint32');
if (byteOrderMark ~= uint32(hex2dec('1A2B3C4D'))) % wrong format
  error('Wrong byte order for this limited script.  FIXME.');
end

%fseek(fid, 8, 'bof'); % skip version and byte-order-mark
numCams = fread(fid, 1, '*uint32');


headers = cell(1,numCams);
cameraNumbers = ones(1,numCams) * (0 / 0);
for i=1:numCams
  pos = 9 + (1024 + 4) * (i-1);
  fseek(fid, pos, 'bof');
  cameraNumbers(i) = fread(fid, 1, '*uint32');
  headers{i} = read_f_header(fid);
end
end

% returns a struct for this header version
function [header] = read_f_header(fid)
version = fread(fid,1, 'int8=>char');
if (version ~= 'f')
  disp(version);
  error('not version "f"')
end

byteOrderMark = fread(fid, 1, '*uint32');
if (byteOrderMark ~= uint32(hex2dec('1A2B3C4D'))) % wrong format
  error('Wrong byte order for this limited script.  FIXME.');
end

frames = fread(fid, 1, '*uint32');
sizeX = fread(fid, 1, '*uint32');
sizeY = fread(fid, 1, '*uint32');
binX = fread(fid, 1, '*uint32');
binY = fread(fid, 1, '*uint32');
acquisitionFrequency = fread(fid, 1, '*uint32');
datetime = fread(fid, 24, 'int8=>char'); % read date string
datetime = transpose(datetime(1:23)); % trim the 0-character and convert to a matlab string
% read null-terminated comment string:
comment = '';
ch = fread(fid, 1, 'int8=>char');
while ch ~= 0
  comment = [comment ch];
  ch = fread(fid, 1, 'int8=>char');
end

frames = double(frames);
sizeX = double(sizeX);
sizeY = double(sizeY);
binX = double(binX);
binY = double(binY);
acquisitionFrequency = double(acquisitionFrequency) * 1e-3;
blockSize = 2*sizeX*sizeY + 8;


header.version = version;
%header.machineformat = machineformat;
header.datetime = datetime;
header.comment = comment;
header.frames = frames;
header.sizeX = sizeX;
header.sizeY = sizeY;
header.binX = binX;
header.binY = binY;
header.acquisitionFrequency = acquisitionFrequency;
header.blockSize = blockSize;
end
