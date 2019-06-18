function [new_data, new_times] = normalizeTimedData(data, times, fps, ...
    varargin)
%NORMALIZETIMEDATA Normalize timed data to a single fps.
% Note: camera clocks are variable but, we need a standard fps to align
% timed data. This function normalizes timed data to a single fps &
% reassigns NaN data to nearest-neighbor time in the downsampled data.
%
%   DATA = NORMALIZETIMEDDATA(DATA, TIME_X, FPS)
%   DATA = NORMALIZETIMEDDATA(DATA, TIME_X, FPS, INTERP_TYPE)
%
%   Inputs:
%   data  - the timed data to normalize (up to 3-dimensional, with time as
%           the last dimension). 
%   times - the timing for the data
%   fps   - the new frames/second for normalization
%   interp_type - the type of interpolation to use
%       default = 'pchip'
%
%   Outputs:
%   data  - the normalized timed data

% What type of interpolation are we using?
interp_type = 'pchip';
if ~isempty(varargin)
    interp_type = varargin{1};
end

% Compute the new time axis.
spf = 1/fps;
new_times = 0:spf:times(end);

% Normalize the data to the new time axis.
switch ndims(data)
    case 1
        new_data = normData(data, times, new_times, interp_type);
    case 2
        new_data = nan(size(data,1),length(new_times));
        for i = 1:size(data,1)
            new_data(i,:) = normData(squeeze(data(i,:)), ...
                times, new_times, interp_type);
        end
    case 3
        new_data = nan(size(data,1),size(data,2),length(new_times));
        for i = 1:size(data,1)
            for j = 1:size(data,2)
                new_data(i,j,:) = normData(squeeze(data(i,j,:)), ...
                    times, new_times, interp_type);
            end
        end
    otherwise
        error('Data is more than 3-dimensional!');
end

end

% Normalize the data to the new time axis.
function new_data = normData(old_data, old_times, new_times, interp_type)

% Is there enough non-NaN data?
isnan_data = isnan(old_data);
if length(old_data) - sum(isnan_data) <= 2
    new_data = nan(1,length(new_times));
    return;
end

% Find gaps of NaNs.
nan_gaps = [];
isnan_gap = false;
for i = 1:length(isnan_data)
    
    % Found the start of a NaN gap.
    if ~isnan_gap && isnan_data(i)
        isnan_gap = true;
        
        % What's the nearest new time?
        % Note: if it's equidistant we bias towards earlier times.
        [~, min_i] = min(abs(new_times - old_times(i)));
        nan_gaps(end+1,1) = min_i;
    end
    
    % Found the end of the NaN gap.
    if isnan_gap && ~isnan_data(i)
        isnan_gap = false;
        
        % What's the nearest new time?
        % Note: if it's equidistant we bias towards later times.
        [~, min_i] = min(abs(flip(new_times) - old_times(i)));
        nan_gaps(end,2) = length(new_times) - min_i + 1;
    end    
end

% Did we end with a NaN?
if isnan_gap
    nan_gaps(end,2) = length(new_times);
end

% Normalize the data to the new time axis.
new_data = interp1(old_times, old_data, new_times, interp_type);

% Fill in the old NaN gaps.
for i = 1:size(nan_gaps,1)
    new_data(nan_gaps(i,1):nan_gaps(i,2)) = nan;
end
end
