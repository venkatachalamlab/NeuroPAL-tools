function [traces, states] = digitizeTraces(traces, varargin)
%DIGITIZETRACES Digitize the neural traces.
%
%   TRACES = DIGITIZETRACES(TRACE)
%   TRACES = DIGITIZETRACES(TRACE, STD_SCALE)
%   TRACES = DIGITIZETRACES(TRACE, STD_SCALE, MIN_FRAMES)
%
%   Inputs:
%   traces     - the neural traces to digitize
%   std_scale  - scalar for digitizing traces by their standard deviation
%                default: 3
%   min_frames - minimum # of frames a state must be held to count as real
%                default: 1
%
%   Outputs:
%   traces     - the digitized neural traces
%   states     - the consecutive frame counts for neural states
%                basal = basal state
%                down = below basal (inhibited)
%                up = above basal (excited)

% The scalar for digitizing traces by their standard deviation.
std_scale = 3;
if ~isempty(varargin)
    std_scale = varargin{1};
end

% The minimum number of frames a state must be held to count as real.
min_frames = 1;
if length(varargin) > 1
    min_frames = varargin{2};
end

% Digitize the traces using the their scaled standard deviation.
std_traces = std_scale * nanstd(traces, 0, 2);
dtraces = floor(traces ./ std_traces);
dtraces = dtraces - mode(dtraces,2);
dtraces(dtraces > 1) = 1;
dtraces((dtraces > -1) & (dtraces < 1)) = 0;
dtraces(dtraces < -1) = -1;

% Count the consecutive frames for neural states & remove spurious states.
states(size(dtraces,1)).down = [];
states(size(dtraces,1)).basal = [];
states(size(dtraces,1)).up = [];
prev_state = NaN;
state = NaN;
count = 0;
for i = 1:size(dtraces,1)
    for j = 1:size(dtraces,2)
        
        % Did we change state?
        if state ~= dtraces(i,j)
            
            % Remove spurious states.
            if j > 1 && count < min_frames
                dtraces(i,(j - count):(j - 1)) = prev_state;
            
            % Record the previous state.
            elseif ~isnan(prev_state) && ~isnan(state) ...
                    && ~isnan(dtraces(i,j))
                switch state
                    case -1
                        states(i).down(end+1) = count;
                    case 0
                        states(i).basal(end+1) = count;
                    case 1
                        states(i).up(end+1) = count;
                end
            end
            
            % Start counting frames for the new state.
            prev_state = state;
            state = dtraces(i,j);
            count = 1;
            
        % Count the number of frames for the state.
        else
            count = count + 1;
        end
    end
end

% Transform the traces to rise (1), hold (0), & fall (-1) states.
traces(:,1) = NaN;
traces(:,2:end) = diff(dtraces,1,2);
traces(traces > 1) = 1;
traces(traces < -1) = -1;
end
