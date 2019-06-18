function peaks = findTracePeaks(traces, fps, varargin)
%FINDTRACEPEAKS Find peak activity in the neural traces.
%
%   Inputs:
%   traces - the neural traces to find peaks in
%   fps    - the sampling rate in frames/second
%
%   Outputs:
%   peaks  - the peaks in the neural traces

std_thresh = 1;
time_thresh = 3;

med_traces = nanmedian(traces, 2);
traces = traces - med_traces;
min_height = std_thresh * nanstd(traces, 0, 2); % + med_traces;
min_dist = round(time_thresh * fps);

peaks = [];
x = (0:(size(traces,2) - 1))/fps;
for i = 1:size(traces,1)
    [~, peaks(i).on] = findpeaks(traces(i,:), ...
        'MinPeakHeight', min_height(i), 'MinPeakDistance', min_dist);
    [~, peaks(i).off] = findpeaks(-traces(i,:), ...
        'MinPeakHeight', min_height(i), 'MinPeakDistance', min_dist);
    figure;
    plot(x,traces(i,:), 'g');
    hold on;
    plot(peaks(i).on/fps, traces(i, peaks(i).on), 'bv');
    plot(peaks(i).off/fps, traces(i, peaks(i).off), 'rv');
    xlabel('Time (s)');
    ylabel('F (a.u.)');
    set(gcf, 'Position', [16, 600, 1900, 300]);
    axis tight;
end
end

