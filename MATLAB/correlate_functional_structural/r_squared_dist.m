function r2 = r_squared_dist(row1, many_rows)
% R_SQUARED coefficient of variation for two vectors, suitable for pdist

% Initialize the data sizes.
N = size(row1, 2);
M2 = size(many_rows, 1);
r2 = zeros(M2, 1);

% Determine the data points.
row1_data = ~isnan(row1);
for i = 1:size(many_rows, 1)
    
    % Comparison row.
    row2 = many_rows(i,:);
    row2_data = ~isnan(row2);
    
    % Compare the rows.
    data = row1_data & row2_data;
    r = corrcoef(row1(data),row2(data));
    r2(i) = 1-r(2)^2;

end
end