function pairs = ranked_pairs_from_pdist(P)
% pairs = RANKED_PAIRS_FROM_PDIST(P)
%
%   Returns an Nx3 array of X-row-index-pairs for the lowest to highest
%   distance values returned by P = pdist(X). NaN values are ignored, so N
%   may be less than length(P);
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

[~,I] = sort(P);

pairs = nan(0, 3);
S = triu(squareform(P));

for i = 1:length(P)
    val = P(I(i));
    if isnan(val)
        continue
    end

    lin_idxs = find(S==val);

    for j = 1:length(lin_idxs)
        lin_idx = lin_idxs(j);
        [row, col] = ind2sub(size(S), lin_idx);
        pairs(end+1,:) = [row, col, val];
    end

end