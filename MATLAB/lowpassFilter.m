function traces = lowpassFilter(traces, low_N, crop_N)
%LOWPASSFILTER Clean the traces with a lowpass filter.
%
%   Input:
%       traces - the traces to filter
%       low_N - the number of lowpass samples to use
%       crop_N - the number of samples to crop
%   Output:
%       traces - the filtered traces

% Build the lowpass filter.
filterTrace = @(x) ...
    lowpass(...
    crop_front(interp_nans(x, 'previous'), crop_N), ...
    low_N);

% Filter the traces.
for i = 1:size(traces,1)
    traces(i,:) = filterTrace(traces(i,:));
end
end
