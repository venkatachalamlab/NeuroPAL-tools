function makeMIPMovie(filename, varargin)
%MAKEMIPMOVIE Make a maximum intensity projection movie.
%
%   MAKEMIPMOVIE(FILENAME)
%   MAKEMIPMOVIE(FILENAME, FPS)
%   MAKEMIPMOVIE(FILENAME, FPS, SCALE)
%   MAKEMIPMOVIE(FILENAME, FPS, SCALE, FRAMES)
%   MAKEMIPMOVIE(FILENAME, FPS, SCALE, FRAMES, SPEED)
%   MAKEMIPMOVIE(FILENAME, FPS, SCALE, FRAMES, SPEED, COMPRESSION)
%
%	Inputs:
%   filename    = the file containing the MIP data for the movie.
%   fps         = the frames/second for the movie
%       default: 4
%   scale       = the scaling factor for the movie height & width
%       default: 4
%   frames      = which frames should we show?
%       default: 1:length(frames) (show every frame -- use [] for default)
%   speed       = the playback speed for the movie
%       default: 1 (1x playback)
%   compression = crude compression rate, implemented by dropping frames
%       default: 1 (keep all frames)

% What is the frames/second for the movie?
fps = 4;
if ~isempty(varargin)
    fps = varargin{1};
end

% Are we scaling the movie height & width.
scale = 4;
if length(varargin) > 1
    scale = varargin{2};
end

% Which frames should we show?
frames = [];
if length(varargin) > 2
    frames = varargin{3};
end

% At what speed should we show the movie?
speed = 1;
if length(varargin) > 3
    speed = varargin{4};
end

% Are we compressing the movie (dropping frames)?
compression = 1;
if length(varargin) > 4
    compression = varargin{5};
    
    % Sanity check the compression.
    compression = round(compression);
    if compression < 1
        compression = 1;
    end
end

% Load the frames/second.
%fps = [];
%fps_filename = strrep(filename, '.mat', '_traces.mat');
%load(fps_filename, 'fps');
%if isempty(fps)
%    error('"%s" does not contain the "fps" for the movie!', fps_filename);
%end

% Load the movie data.
data = [];
load(filename, 'data');
if isempty(data) || ndims(data) ~= 3
    error('The "data" in "%s" cannot be converted to a MIP movie!', ...
        filename);
end

% Which frames are we showing?
if isempty(frames)
    frames = 1:compression:size(data,3);
    
% Check the bounding frames & compression.
else
    
    % Check the bounding frames.
    frames(1) = max(frames(1), 1);
    frames(end) = min(frames(end), size(data,3));
    
    % Check the compression.
    if length(frames) > 1
        compression = round(mean(diff(frames)));
    end
end

% Scale the movie data.
data = data(:,:,frames);
if scale ~= 1
    data = uint8(imresize(data, scale));
end
movie_data(:,:,1,:) = data;

% Make the movie.
movie_filename = strrep(filename, '.mat', '');
movie = VideoWriter(movie_filename, 'MPEG-4');
movie.FrameRate = fps * (speed / compression);
open(movie);
writeVideo(movie, movie_data);
close(movie);
end
