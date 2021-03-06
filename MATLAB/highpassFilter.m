function traces = highpassFilter(traces, high_N, crop_N)
%HIGHPASSFILTER Clean the traces with a highpass filter.
%
%   Input:
%       traces - the traces to filter
%       high_N - the number of highpass samples to use
%       crop_N - the number of samples to crop
%   Output:
%       traces - the filtered traces

% Build the acausal filter.
filterTrace = @(x) ...
    highpass(...
    crop_front(interp_nans(x, 'previous'), crop_N), ...
    high_N);

% Filter the traces.
for i = 1:size(traces,1)
    traces(i,:) = filterTrace(traces(i,:));
end
end
