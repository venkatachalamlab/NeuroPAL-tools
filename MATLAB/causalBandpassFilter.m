function traces = causalBandpassFilter(traces, low_N, high_N, crop_N)
%CAUSALBANDPASSFILTER Clean the traces with an causal bandpass filter.
%
%   Input:
%       traces - the traces to filter
%       low_N - the number of lowpass samples to use
%       high_N - the number of highpass samples to use
%       crop_N - the number of samples to crop
%   Output:
%       traces - the filtered traces

% Build the acausal filter.
filterTrace = @(x) ...
    highpass(...
    lowpass(...
    crop_front(interp_nans(x, 'previous'), crop_N), ...
    low_N), ...
    high_N);

% Filter the traces.
for i = 1:size(traces,1)
    traces(i,:) = filterTrace(traces(i,:));
end
end
