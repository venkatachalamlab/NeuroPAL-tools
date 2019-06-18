%% DON'T EDIT HERE!!!
%% Initilaize user-dependent variables.

% Are we debugging?
is_debug = false;

% Subdirectory for sidedeness of imaging volumes.
left = 'Left Side/';
right = 'Right Side/';
ventral = 'Ventral/';

% Animal segment.
head = 'head';
tail = 'tail';

% Rotate the volume 180 degrees?
isFlip = false;

% Color intensity boost.
mNeptune_boost = 1;
BFP_boost = 1;
CyOFP_boost = 1;
TagRFP_boost = 1;

% Remove background noise?
% Note: weak channels like mNeptune suffer as a result.
isRemovemNeptuneBackground = false;
isRemoveCyOFPBackground = false;
isRemoveBFPBackground = false;
isRemoveTagRFPBackground = false;

%% SET THESE!!!
%% User-dependent variables.

% Near side of the worm = left, right, or ventral.
%side_dir = ventral;
%side_dir = left;
side_dir = right;

% Animal name & number.
animal_number = 41;
image_volume_name = 'id302';
%image_volume_name = 'id_tail_304';

% Animal segment = head or tail.
animal_seg = head;
%animal_seg = tail;

% Color adjustment.
mNeptune_adj = [0, 23200, 3.2];
CyOFP_adj = [0, 4500, 1.2];
BFP_adj = [0, 2400, 1];
TagRFP_adj = [0, 7900, 1.6];

%% DON'T EDIT BELOW THIS LINE!!!
%% Construct user-dependent variables.

% Rotate the volume 180 degrees?
if strcmp(side_dir,left)
    isFlip = true;
end

% Root directory for imaging volumes.
root_dir = '/Users/ev/Desktop/Vivek Worms to ID/';

% Construct the animal input & ouput dirs.
animal_dir = ['animal_' num2str(animal_number,'%03d') '/'];
animal_name = ['animal_' num2str(animal_number,'%03d') '_' animal_seg];

%% Image volume variables.
addpath('./deps');
root = [root_dir side_dir animal_dir];
filename = fullfile(root, sprintf('%s.nd2', image_volume_name));
tiff_filename = fullfile(root, sprintf('%s.ome.tiff', image_volume_name));

asset_directory = fullfile(root, 'annotator_jpegs');

output_directory = fullfile(asset_directory, animal_name);
make_directory(output_directory);

empty_slices = [1:2, 20:21];

%%

original = get_array_data(BioFormatsArray(filename));

% Dimension 4 of A matches nd2 order:
%  1 mNeptune
%  2 RFP emitted into Red
%  3 BFP
%  4 RFP Secondary [empty for older data]
%  5 cyOFP
%  6 RFP Primary

%% Fix orientation

A = original;

if isFlip
    A = flip(A, 1);
    A = flip(A, 3);
end

%% Utilities and filtering parameters

S = size(A);

mip_z = @(img, channel) max_intensity_z(img(:,:,:,channel));

background_filter = fspecial('gaussian', 15, 4);
get_avg_background = @(x) uint16(imfilter(mean(x,3), background_filter));

crop_border = 5; % pixels

cb_ = crop_border - 1;
crop_top = @(x) subsasgn(x, struct(...
    'type', '()', ...
    'subs', {{1:cb_+1, ':', ':', ':'}}), 0);
crop_bottom = @(x) subsasgn(x, struct(...
    'type', '()', ...
    'subs', {{S(1)-cb_:S(1), ':', ':', ':'}}), 0);
crop_left = @(x) subsasgn(x, struct(...
    'type', '()', ...
    'subs', {{':', 1:cb_+1, ':', ':'}}), 0);
crop_right = @(x) subsasgn(x, struct(...
    'type', '()', ...
    'subs', {{':', S(2)-cb_:S(2), ':', ':'}}), 0);

crop = @(x) crop_top(crop_bottom(crop_left(crop_right(x))));


% bandpass filtering
bump_finder = bump_detector([1,1], [4,4], [51 51]);

dc = @(x) deconvlucy(x, fspecial('gaussian', 5, 2), 5);

medflt = @(x, n) medfilt3(x, [5 5 3]);

flt = @(x, h) imfilter(medflt(medflt(dc(double(x)))), h, 'circular');

flt_and_normalize = @(x, bg) autoscale(crop(uint16(flt(x, bump_finder))));

% Make an array like A with four channels: mNep, cyOFP, BFP, all
B = A;
B = crop(B);
B = B(:,:,:,1:4);


%% Prepare mNeptune

mNeptune = A(:,:,:,1);

% Adjust the image.
mNeptune = imadjust3(mNeptune, ...
    [mNeptune_adj(1)/double(65535); mNeptune_adj(2)/double(65535)], ...
    [], mNeptune_adj(3));
mNeptune = uint16(mNeptune);

% Remove the background.
if isRemovemNeptuneBackground
    mNeptune_background = mNeptune(:,:,empty_slices, 1);
    avg_background = get_avg_background(mNeptune_background);
    %new_mNeptune = flt_and_normalize(mNeptune - avg_background);
    new_mNeptune = autoscale(crop(mNeptune - avg_background));
else
    new_mNeptune = mNeptune;
end
    
% Boost the intensity?
if mNeptune_boost ~= 1
    new_mNeptune = new_mNeptune .* mNeptune_boost;
end

% Debugging pics.
if is_debug
    figure(1); clf;
    
    subplot(311);
    imshow(max_intensity_z(mNeptune)*10);
    
    subplot(312);
    imshow(avg_background*10);
    
    subplot(313);
    imshow(max_intensity_z(new_mNeptune));
end

B(:,:,:,1) = new_mNeptune;

%% Prepare cyOFP

CyOFP = A(:,:,:,2);

% Adjust the image.
CyOFP = imadjust3(CyOFP, ...
    [CyOFP_adj(1)/double(65535); CyOFP_adj(2)/double(65535)], ...
    [], CyOFP_adj(3));
CyOFP = uint16(CyOFP);

% Remove the background.
if isRemoveCyOFPBackground
    CyOFP_background = CyOFP(:,:,empty_slices, 1);
    avg_background = get_avg_background(CyOFP_background);
    new_CyOFP = autoscale(crop(CyOFP - avg_background));
else
    new_CyOFP = CyOFP;
end

% Boost the intensity?
if CyOFP_boost ~= 1
    new_CyOFP = new_CyOFP .* CyOFP_boost;
end

% Debugging pics.
if is_debug
    figure(1); clf;
    
    subplot(311);
    imshow(max_intensity_z(CyOFP));
    
    subplot(312);
    imshow(avg_background);
    
    subplot(313);
    imshow(max_intensity_z(new_CyOFP));
end

B(:,:,:,2) = new_CyOFP;

%% Prepare BFP

BFP = A(:,:,:,3);

% Adjust the image.
BFP = imadjust3(BFP, ...
    [BFP_adj(1)/double(65535); BFP_adj(2)/double(65535)], ...
    [], BFP_adj(3));
BFP = uint16(BFP);

% Remove the background.
if isRemoveBFPBackground
    BFP_background = BFP(:,:,empty_slices, 1);
    avg_background = get_avg_background(BFP_background);
    new_BFP = autoscale(crop(BFP - avg_background));
else
    new_BFP = BFP;
end

% Boost the intensity?
if BFP_boost ~= 1
    new_BFP = new_BFP .* BFP_boost;
end

% Debugging pics.
if is_debug
    figure(1); clf;
    
    subplot(311);
    imshow(max_intensity_z(BFP));
    
    subplot(312);
    imshow(avg_background);
    
    subplot(313);
    imshow(max_intensity_z(new_BFP));
end

B(:,:,:,3) = new_BFP;

%% Prepare GCaMP6s

TagRFP = A(:,:,:,4);

% Adjust the image.
TagRFP = imadjust3(TagRFP, ...
    [TagRFP_adj(1)/double(65535); TagRFP_adj(2)/double(65535)], ...
    [], TagRFP_adj(3));
TagRFP = uint16(TagRFP);

% Remove the background.
if isRemoveTagRFPBackground
    TagRFP_background = TagRFP(:,:,empty_slices, 1);
    avg_background = get_avg_background(TagRFP_background);
    new_TagRFP = autoscale(crop(TagRFP - avg_background));
else
    new_TagRFP = TagRFP;
end

% Boost the intensity?
if TagRFP_boost ~= 1
    new_TagRFP = new_TagRFP .* TagRFP_boost;
end

% Debugging pics.
if is_debug
    figure(1); clf;
    
    subplot(311);
    imshow(max_intensity_z(TagRFP));
    
    subplot(312);
    imshow(avg_background);
    
    subplot(313);
    imshow(max_intensity_z(new_TagRFP));
end

B(:,:,:,4) = new_TagRFP;

%% Write Tiff

%SB = size(B);
%K = reshape(B, [SB(1), SB(2), SB(3)*SB(4)]);

%save_tiff_stack(K, tiff_filename);
bfsave(B, tiff_filename);

%%

channels{1} = 'LM';

         %mNep   %cyOFP  %BFP     %PanNeuronal
LUT{1} = [1.0     0.0     0.0      0.0;
          0.0     0.4     0.0      0.0;
          0.0     0.0     1.0      0.0];
     
channels{2} = 'LM+PN';

         %mNep   %cyOFP  %BFP   %RFP1
LUT{2} = [1.0     0.0     0.0      0.6;
          0.0     0.4     0.0      0.6;
          0.0     0.0     1.0      0.6];

channels{3} = 'LM+PN Corrected';

LM = 0.3; %double(1.0)/double(3.0);
PN = 0.25;

         %mNep   %cyOFP  %BFP   %RFP1
LUT{3} = [LM   0.0  0.0   PN;
          0.0  LM   0.0   PN;
          0.0  0.0  LM    PN];
      
clear colors;
colors{1} = get_slice(B,1);
colors{2} = get_slice(B,2);
colors{3} = get_slice(B,3);
colors{4} = get_slice(B,4);

%C = 3.0*get_array_data(ColorArray(colors, 4, LUT{3}, true, 5));

% Add the panneuronal to the color channels
C=B(:,:,:,1:3);
LM_scale = [0.2, 0.8, 0.7]; % [mNep cyOFP BFP]
PN_scale = 1;%2;
for i=1:size(C,3)
    
    % Scale the channels according to background noise.
    D = squeeze(C(:,:,i,:));
    for j = 1:3
        D(:,:,j) = D(:,:,j) * LM_scale(j);
    end
    
    % Lower the panneuronal if it obscures neural coloring.
    max_LM = max(D, [], 3);
    PN_intensity = squeeze(B(:,:,i,4));
    PN_add = PN_intensity - max_LM;
    PN_add = PN_add./PN_scale;
    
    % Scale the panneuronal so it doesn't obscure neural coloring.
    C(:,:,i,:) = C(:,:,i,:) + PN_add;
    %C(:,:,i,:) = C(:,:,i,:) + ((B(:,:,i,4) - (max(C(:,:,i,:), [], 4))./PN_scale));
    %C(:,:,i,:) = C(:,:,i,:) + ((B(:,:,i,4) - uint16(sum(C(:,:,i,:), 4))));
end

figure; imshow(max_intensity_z(C));

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

% Show the TIFF.
system(['/usr/local/bin/ijview ' tiff_filename '&'], '-echo');