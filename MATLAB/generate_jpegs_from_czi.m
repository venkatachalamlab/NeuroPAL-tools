root = 'V:\data\170301 NeuroPAL Volumes for Vivek';
path = 'B';
name = '1_YApSP';
filename = fullfile(root, path, sprintf('%s.czi', name));

asset_directory = fullfile(root, 'generated');

output_directory = fullfile(asset_directory, [path '_' name]);
make_directory(output_directory);

%%
bio_formats_data = bfopen(filename);

metadata = bio_formats_data{2};
omeMeta = bio_formats_data{4};

size_X = omeMeta.getPixelsSizeX(0).getValue; % 512
size_Y = omeMeta.getPixelsSizeY(0).getValue; % 512
size_Z = omeMeta.getPixelsSizeZ(0).getValue; % 31
size_T = omeMeta.getPixelsSizeT(0).getValue; % 615
size_C = omeMeta.getPixelsSizeC(0).getValue; % 2

colors = {};
z_ratio = 1;

for c = 1:size_C
    
    idx = 1;
    
    for i = c:size_C:(size_Z*size_C)
        
        for j = 1:z_ratio
            colors{c}(:,:,idx) = bio_formats_data{1}{i,1};
            idx = idx + 1;
        end
        
    end
    
end

%%

       %PN     %DIC    %BFP    %mNep   %CyOFP
LUT = [0.5     0.9     0.0     0.0     1.0;
       0.5     0.9     0.0     1.0     0.0;
       0.5     0.9     1.0     0.0     0.0];

C = get_array_data(ColorArray(colors, 4, LUT, true, 5));

figure(1); clf; imshow(max_intensity_z(C));


%%

make_filename = @(view, view_idx, channel_idx, t_idx) fullfile(...
    output_directory, ...
    sprintf('%s_%d_%d_%d.jpg', view, view_idx, channel_idx, t_idx));

write = @(x, view, view_idx, channel_idx, t_idx) imwrite(uint8(x/2^8), ...
    make_filename(view, view_idx, channel_idx, t_idx), ...
    'Quality', 50);

write_x = @(im, idx) write(im, 'X', idx, 0, 0);
write_y = @(im, idx) write(im, 'Y', idx, 0, 0);
write_z = @(im, idx) write(im, 'Z', idx, 0, 0);

% Max intensity projections

MIPZ = max_intensity_z(C);
MIPX = max_intensity_x(C);
MIPY = max_intensity_y(C);

write(MIPX, 'MIPX', 0, 0, 0)
write(MIPY, 'MIPY', 0, 0, 0)
write(MIPZ, 'MIPZ', 0, 0, 0)

% Slices

% Z
for i = 1:size(C, 3)
    
    slice = squeeze(C(:,:,i,:));
    write_z(slice, i);
    
end

% Y
for i = 1:size(C, 1)
    
    slice = squeeze(C(i,:,:,:));
    slice = permute(slice, [2 1 3]);
    write_y(slice, i);
    
end

% X
for i = 1:size(C, 2)
    
    slice = squeeze(C(:,i,:,:));
    write_x(slice, i);
    
end

