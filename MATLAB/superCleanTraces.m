function [traces, removed] = superCleanTraces(traces, fps, varargin)
%UNTITLED2 Summary of this function goes here
%
%   Detailed explanation goes here

% Remove traces with too little data.
nans_threshold = 0.1;
num_data = sum(~isnan(traces),2);
remove_nans = num_data < size(traces,2) * nans_threshold;
traces(remove_nans,:) = nan;

% Filter out high-frequency noise.
% Butterworth (flat passband), 2nd order, & 0-phase.
interp_method = 'pchip'; %'linear';
filt_order = 70; %12; %18; % 3rd order
half_freq = 0.2;
cut_freq = 1/3; % 1/3Hz cutoff.
norm_freq = cut_freq / (fps/2);
%df = designfilt('lowpassiir', 'FilterOrder', filt_order, ...
%    'HalfPowerFrequency', half_freq, 'SampleRate', fps, ...
%    'DesignMethod', 'butter');
%df = designfilt('lowpassfir', 'FilterOrder', filt_order, ...
%    'CutoffFrequency', norm_freq);
%df = designfilt('lowpassfir','PassbandFrequency',1/8,...
%  'StopbandFrequency',1/4,'FilterOrder',80, ...
%  'SampleRate', fps, 'DesignMethod', 'ls');
df = designfilt('lowpassfir','PassbandFrequency',1/4,...
  'StopbandFrequency',1/3,'FilterOrder',12, ...
  'SampleRate', fps, 'DesignMethod', 'ls');
keepI = find(~remove_nans);
for i = 1:length(keepI)
    y = traces(keepI(i),:);
    y_data = ~isnan(y);
    y = fillmissing(y, interp_method);
    y = filtfilt(df,y);
    traces(keepI(i),y_data) = y(y_data);
end

% Remove weak traces.
min_threshold = 16;
max_traces = nanmax(traces, [], 2);
remove_mins = max_traces < min_threshold;
traces(remove_mins,:) = nan;

% Remove traces with weak signal.
std_threshold = 4;
std_traces = nanstd(traces, 0, 2);
remove_stds = std_traces < std_threshold;
traces(remove_stds,:) = nan;

% Clean up the traces.
sgolay_scale = 1;
gcamp_signal_tau = 3;
outlier_sigma = 6;
sgolay_win = sgolay_scale * (gcamp_signal_tau * fps);
traces = cleanTraces(traces, outlier_sigma, 2, [], 'sgolay', sgolay_win);
traces = cleanTraces(traces, outlier_sigma, 2);

% Which neurons were removed?
removed = remove_nans | remove_mins | remove_stds;
end
