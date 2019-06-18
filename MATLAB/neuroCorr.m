function [cor, lag, neuroI] = neuroCorr(traces, frames)
%NEUROCORR Compute unbiased neural cross-correlations, pairwise.
%
%   [COR, LAG, NEUROI] = NEUROCORR(TRACES, FRAMES)
%
%   Inputs:
%   traces - the traces for computing correlations
%   frames - the frames within the traces for computing correlations
%
%   Outputs:
%   cor    - the pairwise, neural cross-correlations
%   lag    - the cross-correlation lags
%   neuroI - the indices for neuron pairs, per cross-correlation

% Preallocate memory.
num_neurons = size(traces, 1);
num_pairs = nchoosek(num_neurons, 2);
num_corrs = 2 * length(frames) - 1;
neuroI = nan(num_pairs, 2);
cor = nan(num_pairs, num_corrs);
lag = nan(num_pairs, num_corrs);

% Compute the windowed correlations.
neuroI_ij = 1;
for i = 1:num_neurons
    for j = i:num_neurons
        
        % Organize the neural traces.
        n1 = traces(i, frames);
        n2 = traces(j, frames);
        
        % Compute the pairwise correlation.
        neuroI(neuroI_ij,:) = [i,j];
        [cor(neuroI_ij,:), lag(neuroI_ij,:)] = neuroCorr1(n1, n2);
        
        % Increment the count.
        neuroI_ij = neuroI_ij + 1;
    end
end
end

%% Compute the pairwise neural corraltions.
function [cor, lag] = neuroCorr1(n1, n2)

% Preallocate memory.
num_corrs = 2 * length(n1) - 1;
cor = nan(1, num_corrs);
lag = nan(1, num_corrs);

% Compute the positive-lag correlations.
for i = 0:(length(n1)-1)
    
    % Compute the lag.
    corrI = length(n1) + i;
    lag(corrI) = i;

    % Compute the positive correlation.
    corrs = n1((i + 1):end) .* n2(1:(end - i));
    cor(corrI) = nansum(corrs) / sum(~isnan(corrs));
end

% Compute the negative-lag correlations.
for i = 1:(length(n1)-1)
    
    % Compute the lag.
    corrI = length(n1) - i;
    lag(corrI) = -i;

    % Compute the positive correlation.
    corrs = n2((i + 1):end) .* n1(1:(end - i));
    cor(corrI) = nansum(corrs) / sum(~isnan(corrs));
end
end
