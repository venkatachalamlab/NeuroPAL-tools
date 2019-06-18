function index = nanMedianIndex(data, dim)
%NANMEDIANVECTOR Find the index of the median value(s) in an matrix.
%
%   INDEX = NANMEDIANINDEX(DATA)
%   INDEX = NANMEDIANINDEX(DATA, DIM)
%
%   Inputs:
%   data = the data
%   dim  = the dimension in which to find the median value(s)
%          default: 1
%
%   Outputs:
%   index = the index of the median value(s)

med = nanmedian(data, dim);
[~, index] = min(nanmean(abs(data - med)));
end
