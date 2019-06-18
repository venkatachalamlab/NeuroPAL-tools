%% Combine MIPs with color movie.

% Initialize the color data & GCaMP MIP.
color_data_file = "id101.ome.tiff";
mip_data_file = "run101_mip.mat";
time_data_file = "run101times.mat";

% Determine the frames/s.
TM = load(time_data_file);
fps = 1/mean(diff(TM.times.times));

% Load the color ID image.
color_data = load_tiff_stack(color_data_file);
S = size(color_data);
color_data = reshape(color_data, [S(1), S(2), 21, 4]);

R = color_data(:,:,:,1);
G = color_data(:,:,:,2);
B = color_data(:,:,:,3);

% Clean up haze in red channel (high pass filter)
R64 = double(R)/2^16;
R64 = R64 - smooth3(R64, 'gaussian', [31 31 11], 15);
R = uint16(2^16*R64);
% done.

threshold = @(x, p) x - quantile(x(:), p);
median_filter = @(x, Nx, Nz) medfilt3(x, [Nx, Nx, Nz]);
smooth = @(x) smooth3(double(x)/double(quantile(x(:), 0.9999)), 'gaussian');
process = @(x) uint8(255*smooth(median_filter(threshold(x, 0.7), 5, 3)));
mipz = @max_intensity_z;

combine = @(R, G, B) cat(3, mipz(R), mipz(G), mipz(B));
show3D = @(x) imshow(mipz(x));
show3DRGB = @(R, G, B) imshow(cat(3, mipz(R), mipz(G), mipz(B)));

colorscale = [1.5, 1.0, 1.5];

R_out = process(R) * colorscale(1);
G_out = process(G) * colorscale(1);
B_out = process(B) * colorscale(1);

figure(6); clf; show3DRGB(R_out, G_out, B_out);

% Upsample the image, it's too small.
resize_scale = 2;
R_out = imresize(R_out,resize_scale);
G_out = imresize(G_out,resize_scale);
B_out = imresize(B_out,resize_scale);

% Load the GCaMP MIP frames.
M = load(mip_data_file);
mipg = M.data;

% Smooth the GCaMP MIP frames.
%mipg = mipg - median(mipg,3);
nhood = [4,4];
for i = 1:size(mipg,3)
    mipg(:,:,i) = medfilt2(mipg(:,:,i), nhood);
end

% Translate the GCaMP MIP frames.
start_trans = [3, -1];
end_trans = [4, -1];
trans(1,:) = linspace(start_trans(1), end_trans(1), size(mipg,3));
trans(2,:) = linspace(start_trans(2), end_trans(2), size(mipg,3));
for i = 1:size(mipg,3)
    mipg(:,:,i) = imtranslate(mipg(:,:,i), [trans(1,i), trans(2,i)]);
end

% Upsample the GCaMP MIP frames.
mipg = imresize(mipg,resize_scale * 2);
imRGBW = @(f) cat(3, mipz(R_out) + mipg(:,:,f), mipz(G_out) + mipg(:,:,f), mipz(B_out) + mipg(:,:,f));
showRGBW = @(f) imshow(imRGBW(f));

%% Make the video.

% Open the video.
offset = 0;
speed = 12;
filename = strrep(mip_data_file, '_mip.mat', ...
    ['_' num2str(offset) 'off_' num2str(resize_scale) 'scale_' num2str(speed) 'x.mp4']);
writer = VideoWriter(filename,'MPEG-4');
writer.Quality = 100;
writer.FrameRate = fps * speed;
open(writer);

% Offset the movie.
mov = mipg;
if offset > 1
    mov = move(:,:,offset:end);
end
imRGBW = @(f) cat(3, mipz(R_out) + mov(:,:,f), mipz(G_out) + mov(:,:,f), mipz(B_out) + mov(:,:,f));

% Write the frames.
wb = waitbar(0, "Creating movie.");
for t = 1:size(mov,3)
    im = imRGBW(t);
    writeVideo(writer, im);
    waitbar(t/size(mov,3), wb, "Creating movie.");
end
close(writer);

