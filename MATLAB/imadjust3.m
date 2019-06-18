function J = imadjust3(I, varargin)

% Initialize the adjusted image.
J = zeros(size(I));

% Auto adjust.
if isempty(varargin)
    for i = 1:size(I,3)
        %J(:,:,i) = adapthisteq(I(:,:,i), 'Distribution', 'rayleigh');
        J(:,:,i) = imadjust(I(:,:,i));
    end
    %J = I;
    
% Exact adjustment.
elseif length(varargin) == 1
    for i = 1:size(I,3)
        J(:,:,i) = imadjust(I(:,:,i), varargin{1});
    end
elseif length(varargin) == 2
    for i = 1:size(I,3)
        J(:,:,i) = imadjust(I(:,:,i), varargin{1}, varargin{2});
    end
else
    for i = 1:size(I,3)
        J(:,:,i) = imadjust(I(:,:,i), varargin{1}, varargin{2}, varargin{3});
    end
end
end