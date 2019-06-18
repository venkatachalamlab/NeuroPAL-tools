function movie_from_colors(r, g, b, filename, varargin)
% movie_from_colors(r, g, b, filename)
%
% Generates a movie from r, g, and b matrices (each M x N x T).

if exist(filename, 'file')
    delete(filename);
end

if isempty(b)
    if isempty(g)
        b = zeros(size(r), element_class(r));
    else
        b = zeros(size(g), element_class(g));
    end
end

if isempty(r)
    if isempty(b)
        r = zeros(size(g), element_class(g));
    else
        r = zeros(size(b), element_class(b));
    end
end

if isempty(g)
    if isempty(r)
        g = zeros(size(b), element_class(b));
    else
        g = zeros(size(r), element_class(r));
    end
end

r = double(r);
g = double(g);
b = double(b);

writer = VideoWriter(filename,'MPEG-4');
writer.Quality = 100;
writer.FrameRate = 15;
open(writer);

r_max = quantile(reshape(r, [numel(r) 1]), 1-1e-4);
g_max = quantile(reshape(g, [numel(g) 1]), 1-1e-4);
b_max = quantile(reshape(b, [numel(b) 1]), 1-1e-4);

r = r / (r_max+eps);
g = g / (g_max+eps);
b = b / (b_max+eps);

r(r>1) = 1;
g(g>1) = 1;
b(b>1) = 1;

wb = waitbar(0, "Creating movie.");
for t = 1:size(r,3)
    im_r = kron(r(:,:,t), ones(2,2));
    im_g = kron(g(:,:,t), ones(2,2));
    im_b = kron(b(:,:,t), ones(2,2));
    
    im_r = imfilter(im_r, ones(2,2)/4);
    im_g = imfilter(im_g, ones(2,2)/4);
    im_b = imfilter(im_b, ones(2,2)/4);
    
    im = cat(3, im_r, im_g, im_b);
    writeVideo(writer, im);
    waitbar(t/size(r,3), wb, "Creating movie.");
end