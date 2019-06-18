function plotExpMovie(exp_name, varargin)
%PLOTEXPMOVIE Plot an experiment's neural activity as a movie.
%
%   PLOTEXPMOVIE(EXP_NAME)
%   PLOTEXPMOVIE(EXP_NAME, PLOT_MODE)
%   PLOTEXPMOVIE(EXP_NAME, PLOT_MODE, SIDE_MODE)
%   PLOTEXPMOVIE(EXP_NAME, PLOT_MODE, SIDE_MODE, FRAMES)
%   PLOTEXPMOVIE(EXP_NAME, PLOT_MODE, SIDE_MODE, FRAMES, MOVIE_NAME, SPEED)
%   PLOTEXPMOVIE(EXP_NAME, PLOT_MODE, SIDE_MODE, FRAMES, MOVIE_NAME, SPEED,
%                COMPRESSION)
%
%   Inputs:
%   exp_name    = the name of the experiment to plot
%   plot_mode   = how should we plot the neural traces?
%       default: 1
%       1 = plot the neural activity via changes in size & transparency
%       2 = plot the neural activity via changes in size
%       3 = plot the neural activity via changes in transparency
%       4 = don't plot the neural activity
%   side_mode   = which side of the worm should we show?
%       default: 1
%       1 = show both sides of the worm
%       2 = show the left side of the worm
%       3 = show the right side of the worm
%   frames      = which frames should we show?
%       default: 1:length(frames) (show every frame -- use [] for default)
%   movie_name  = save the movie to a file with the name 'movie_name'
%       default: [] (don't save the movie to a file)
%   speed       = the playback speed for the movie
%       default: 1 (1x playback)
%   compression = crude compression rate, implemented by dropping frames
%       default: 1 (keep all frames)

% Are we combining left/right neural designations?
plot_mode = 1;
if ~isempty(varargin)
    plot_mode = varargin{1};
end
switch plot_mode
    case 1
        plot_mode_label = 'Neuron Size & Transparency = Activity';
    case 2
        plot_mode_label = 'Neuron Size = Activity';
    case 3
        plot_mode_label = 'Neuron Transparency = Activity';
    case 4
        % do nothing
    otherwise
        error('Unknown plot mode: %d', plot_mode);
end

% Which side of the worm should we show?
side_mode = 1;
if length(varargin) > 1
    side_mode = varargin{2};
end

% Which frames should we show?
frames = [];
if length(varargin) > 2
    frames = varargin{3};
end

% Should we save the movie to a file and, if so, what is the name?
movie_name = [];
if length(varargin) > 3
    movie_name = varargin{4};
    
    % Does the movie already exist?
    if ~isempty(movie_name) && exist(movie_name, 'file')
        error('"%" already exists!', movie_name);
    end
end

% At what speed should we show the movie?
speed = 1;
if length(varargin) > 4
    speed = varargin{5};
end

% Are we compressing the movie (dropping frames)?
compression = 1;
if length(varargin) > 5
    compression = varargin{6};
    
    % Sanity check the compression.
    compression = round(compression);
    if compression < 1
        compression = 1;
    end
end

%% Initialize orientation, neurons,and stimulus properties.
% Initialize the neuron colors.
load('neuron_colors.mat');
neuron_colors = neurons;
neuron_color_list = {neuron_colors.name};
white = ones(1,3);
panneuronal = white * 0.5;
gray = white * 0.7;
brighten_scale = 1.1;

% Initialize the neuron transparency.
default_neuron_alpha = 0.5;
gcamp_alpha_scale = 0.7;
gcamp_alpha_off = 0.05;

% Initialize the neuron size and shape.
default_neuron_size_scale = 1.5;
neuron_size = 30;
gcamp_size_scale = 2;
gcamp_size_off = 0.5;
[x, y, z] = sphere(neuron_size);
neuron_patch = surf2patch(x,y,z);

% Initialize the padding for plot limits.
pos_lim_pad = 5;

% Initialize the stimuli.
load('stimuli.mat');
stim_list = stimulus;
stim_list_keywords = {stim_list.keyword};

% Initialize the ganglia.
load('ganglia_LR.mat');
ganglion = ganglion_LR;

%% Organize the experiment.
% Open the experiment.
exp_data = open(exp_name);
fps = exp_data.fps;

% Organize the stimuli.
stim = [];
exp_stim = exp_data.stimulus;
for i=1:length(exp_stim)
    
    % Which of the stimuli is being applied?
    stim_name = lower(exp_stim{i,1});
    iStim = find(cellfun(@(x) contains(stim_name,x), stim_list_keywords));
    if isempty(iStim)
        error('Unknown stimulus: %s', exp_stim{i,1});
    end
    
    % Organize the stimulus.
    stim(i).num = iStim;
    stim(i).on = round(exp_stim{i,2} * fps);
    stim(i).off = round(exp_stim{i,3} * fps);
end

% Organize the neuron names.
neuron_names = exp_data.neuron_names;
if iscategorical(neuron_names)
    neuron_names = cellstr(neuron_names);
end
    
% Organize the neuron colors.
neuron_RGBs = zeros(length(neuron_names), 3);
for i=1:length(neuron_names)
    
    % Get the neuron's color.
    neuron_RGBs(i,:) = panneuronal;
    name = stripNeuronLR(neuron_names{i}, true);
    iColor = find(strcmp(neuron_color_list, name));
    if ~isempty(iColor)
        if ~all(neuron_colors(iColor).RGB == 0) % ignore panneuronals
            neuron_RGBs(i,:) = neuron_colors(iColor).RGB;
        end
        
        % Brighten the colors.
        neuron_RGBs(i,:) = min(neuron_RGBs(i,:) * brighten_scale, white);
    end
end
    
% Organize and scale the neural positions.
pos = exp_data.positions;
pos(:,1,:) = pos(:,1,:) .* 0.5;
pos(:,2,:) = pos(:,2,:) .* 0.5;
pos(:,3,:) = pos(:,3,:) .* 1.5;

% Sanity check the positions.
% Note: remove this when the bug is fixed.
max_pos_threshold = 10000;
min_pos_threshold = -10000;
[iNeuron, ~, iFrame] = ind2sub(size(pos), find(pos > max_pos_threshold));
for i=1:length(iNeuron)
    pos(iNeuron(i),:,iFrame(i)) = nan;
end
[iNeuron, ~, iFrame] = ind2sub(size(pos), find(pos < min_pos_threshold));
for i=1:length(iNeuron)
    pos(iNeuron(i),:,iFrame(i)) = nan;
end

% Clean up the neural traces.
gcamp = exp_data.gcamp;
gcamp = cleanTraces(gcamp);
[iNeuron, iFrame] = find(isnan(gcamp));
for i=1:length(iNeuron)
    pos(iNeuron(i),:,iFrame(i)) = nan;
end

% Warn the user about missing frames.
missing_frames = find(all(isnan(gcamp),1));
if ~isempty(missing_frames)
    missing_str = ['Missing frames: ' num2str(missing_frames(1))];
    for i=2:length(missing_frames)
        missing_str = [missing_str ', ' num2str(missing_frames(i))];
    end
    warning(missing_str);
end

% Warn the user about missing neurons.
missing_neurons = find(all(isnan(gcamp),2));
if ~isempty(missing_neurons)
    missing_str = ['Missing neurons: ' neuron_names{missing_neurons(1)}];
    for i=2:length(missing_neurons)
        missing_str = [missing_str ', ' neuron_names{missing_neurons(i)}];
    end
    warning(missing_str);
end

% Compute the neurons' size & transparency.
if plot_mode == 1 || plot_mode == 2
    neuron_sizes = (gcamp * gcamp_size_scale) + gcamp_size_off;
else
    neuron_sizes = ones(size(gcamp)) * default_neuron_size_scale;
end
if plot_mode == 1 || plot_mode == 3
    neuron_alphas = (gcamp * gcamp_alpha_scale) + gcamp_alpha_off;
else
    neuron_alphas = ones(size(gcamp)) * default_neuron_alpha;
end

%% Orient the worm.
% Initialize the head neurons.
ant_head_L = {'OLQDL','OLQVL'};
ant_head_R = {'OLQDR','OLQVR'};
post_head_L = {'AVAL','AVEL'};
post_head_R = {'AVAR','AVER'};
[~, ~, iAHL] = intersect(ant_head_L, neuron_names);
[~, ~, iAHR] = intersect(ant_head_R, neuron_names);
[~, ~, iPHL] = intersect(post_head_L, neuron_names);
[~, ~, iPHR] = intersect(post_head_R, neuron_names);

% Initialize the tail neurons.
ant_tail_L = {'LUAL','PVCL'};
ant_tail_R = {'LUAR','PVCR'};
post_tail_L = {'PVNL','PLML'};
post_tail_R = {'PVNR','PLMR'};
[~, ~, iATL] = intersect(ant_tail_L, neuron_names);
[~, ~, iATR] = intersect(ant_tail_R, neuron_names);
[~, ~, iPTL] = intersect(post_tail_L, neuron_names);
[~, ~, iPTR] = intersect(post_tail_R, neuron_names);

% Are we showing the head?
isHead = false;
if length(iAHL) == length(ant_head_L) && ...
        length(iAHR) == length(ant_head_R) && ...
        length(iPHL) == length(post_head_L) && ...
        length(iPHR) == length(post_head_R)
    isHead = true;
    body_text = 'Head';
    
    % Compute the head neurons' positons.
    %pos_AHL = squeeze(nanmean(pos(iAHL,:,:), 1));
    %pos_AHR = squeeze(nanmean(pos(iAHR,:,:), 1));
    pos_PHL = squeeze(nanmean(pos(iPHL,:,:), 1));
    pos_PHR = squeeze(nanmean(pos(iPHR,:,:), 1));
    
    % Compute the vectors for the worm's orientation.
    %vAPL = pos_PHL - pos_AHL; % AP vector on left
    %vAPR = pos_PHR - pos_AHR; % AP vector on right
    %vAP = (vAPL + vAPR) / 2; % AP vector
    vLR = pos_PHR - pos_PHL; % LR vector
    
% Are we showing the tail?
elseif length(iATL) == length(ant_tail_L) && ...
        length(iATR) == length(ant_tail_R) && ...
        length(iPTL) == length(post_tail_L) && ...
        length(iPTR) == length(post_tail_R)
    body_text = 'Tail';

    % Compute the tail neurons' positons.
    pos_ATL = squeeze(nanmean(pos(iATL,:,:), 1));
    pos_ATR = squeeze(nanmean(pos(iATR,:,:), 1));
    %pos_PTL = squeeze(nanmean(pos(iPTL,:,:), 1));
    %pos_PTR = squeeze(nanmean(pos(iPTR,:,:), 1));
    
    % Compute the vectors for the worm's orientation.
    %vAPL = pos_PTL - pos_ATL; % AP vector on left
    %vAPR = pos_PTR - pos_ATR; % AP vector on right
    %vAP = (vAPL + vAPR) / 2; % AP vector
    vLR = pos_ATR - pos_ATL; % LR vector
        
% Nothing to show.
else
    error('"%s" does not contain enough neurons to show the head or tail!', ...
        exp_name);
end

% Rotate the worm so dorsal-ventral = +/- z-axis.
% The microfluidic orients the worm correctly in AP.
% Assume the worm is oriented so AP = x-axis and LR should = y-axis.
% Compute the x-axis rotation using the angle between LR & the y-axis.
x = [1,0,0];
y = [0,-1,0];
%vLR_med = nanmedian(vLR,2); % median LR vector
vLR = vLR ./ sqrt(nansum(vLR.^2,1)); % LR unit vector
vLR_med = vLR(:,nanMedianIndex(vLR,2)); % median LR vector
cross_vLR_y = cross(vLR_med, y);
x_rot_sign = sign(dot(cross_vLR_y, x));
x_rot = rotx(x_rot_sign * atan2d(norm(cross_vLR_y), dot(vLR_med,y)));
for i=1:size(pos,3)
    pos(:,:,i) = (x_rot*pos(:,:,i)')';
end

% Re-compute the head vectors from the rotated worm.
if isHead
    
    % Compute the head neurons' positons.
    pos_AHL = squeeze(nanmean(pos(iAHL,:,:), 1));
    pos_AHR = squeeze(nanmean(pos(iAHR,:,:), 1));
    pos_PHL = squeeze(nanmean(pos(iPHL,:,:), 1));
    pos_PHR = squeeze(nanmean(pos(iPHR,:,:), 1));
    
    % Compute the vectors for the worm's orientation.
    vAPL = pos_PHL - pos_AHL; % AP vector on left
    vAPR = pos_PHR - pos_AHR; % AP vector on right
    vAP = (vAPL + vAPR) / 2; % AP vector
    vLR = pos_PHR - pos_PHL; % LR vector
    
% Re-compute the tail vectors from the rotated worm.
else
    
    % Compute the tail neurons' positons.
    pos_ATL = squeeze(nanmean(pos(iATL,:,:), 1));
    pos_ATR = squeeze(nanmean(pos(iATR,:,:), 1));
    pos_PTL = squeeze(nanmean(pos(iPTL,:,:), 1));
    pos_PTR = squeeze(nanmean(pos(iPTR,:,:), 1));
    
    % Compute the vectors for the worm's orientation.
    vAPL = pos_PTL - pos_ATL; % AP vector on left
    vAPR = pos_PTR - pos_ATR; % AP vector on right
    vAP = (vAPL + vAPR) / 2; % AP vector
    vLR = pos_ATR - pos_ATL; % LR vector
end

% Normalize to unit vectors.
vAP = vAP ./ sqrt(nansum(vAP.^2,1)); % AP unit vector
vLR = vLR ./ sqrt(nansum(vLR.^2,1)); % LR unit vector

% Orient LR orthogonal to AP by removing its AP projection.
vLR = vLR - dot(vLR, vAP) .* vAP;
vLR = vLR ./ sqrt(nansum(vLR.^2, 1));

% Compute the median AP and LR vectors.
%vAP_med = nanmedian(vAP,2);
%vLR_med = nanmedian(vLR,2);
%vAP_med = vAP_med ./ sqrt(nansum(vAP_med.^2,1));
%vLR_med = vLR_med ./ sqrt(nansum(vLR_med.^2,1));
vAP_med = vAP(:,nanMedianIndex(vAP,2));
vLR_med = vLR(:,nanMedianIndex(vLR,2));

% Orient DV orthogonal to AP and LR.
vDV_med = cross(vLR_med, vAP_med);

% Correct the median vectors to be the x,y,z axes.
rot = [vAP_med'; vLR_med'; vDV_med'];
vAP = rot*vAP;
vLR = rot*vLR;
vDV = -cross(vLR, vAP); % mirror the z-axis to match image conventions

% Re-compute the median vectors for the rotated axes.
%vAP_med = nanmedian(vAP,2);
%vLR_med = nanmedian(vLR,2);
%vDV_med = nanmedian(vDV,2);
%vAP_med = vAP_med ./ sqrt(nansum(vAP_med.^2,1));
%vLR_med = vLR_med ./ sqrt(nansum(vLR_med.^2,1));
%vDV_med = vDV_med ./ sqrt(nansum(vDV_med.^2,1));
vAP_med = vAP(:,nanMedianIndex(vAP,2));
vLR_med = vLR(:,nanMedianIndex(vLR,2));

% Orient LR orthogonal to AP by removing its AP projection.
vLR_med = vLR_med - dot(vLR_med, vAP_med) .* vAP_med;
vLR_med = vLR_med ./ sqrt(nansum(vLR_med.^2, 1));

% Orient DV orthogonal to AP and LR.
vDV_med = -cross(vLR_med, vAP_med); % mirror the z-axis to match image conventions

% Orient the worm using a 3D PCA: rotating to the neural eigen basis.
% Initially this seemed like a good idea. The majority of the neurons
% are distributed anterior-posterior in X, forming the 1st component.
% Thereafter, the neurons are mostly distributed left-right, forming
% the 2nd orthogonal component. And, lastly, the neurons are
% distributed dorsal-ventral, forming the 3rd orthogonal component.
% But, as life would have it, the worm's left and right sides are not
% symmetric. Ganglia on either side are displaced, in X, relative to
% each other. Moreover, the PAG is usually an odd asymmetric shape.
% Consequently, the left-right gangliar displacement often tilts the
% tail towards one side and the irregular shape of the PAG further
% compounds the problem in various axes.
%
% Here's what we tried.
%
% Use every visible neuron or a subset of tail neurons.
% [~, ~, iNeurons] = intersect(tailNeurons, neuron_names, 'stable');
% Compute the worm orientation.
% for i=1:size(pos,3)
%    comps = pca(pos(:,:,i));
%    comps = pca(pos(iNeurons,:,i));
%    vAP(:,i) = comps(:,1);
%    vLR(:,i) = comps(:,2);
% end
%
% Compute the median axes.
% vAP_med = median(vAP, 2);
% vLR_med = median(vLR, 2);
%
% Correct the vector basis to be orthogonal.
% vAP_med = vAP_med ./ sqrt(sum(vAP_med.^2, 1));
% vLR_med = vLR_med - dot(vLR_med, vAP_med) * vAP_med;
% vLR_med = vLR_med ./ sqrt(sum(vLR_med.^2, 1));
% vDV_med = cross(vLR_med, vAP_med);
%
% Orient the tail axes correctly.
%aAP = atan2(norm(cross(vAP_med,vAP_pca_med)), dot(vAP_med,vAP_pca_med));
%aLR = atan2(norm(cross(vLR_med,vLR_pca_med)), dot(vLR_med,vLR_pca_med));
%if aAP < -pi/2 || aAP > pi/2
%    vAP = -vAP;
%    vDV = -vDV;
%end
%if aLR < -pi/2 || aLR > pi/2
%    vLR = -vLR;
%    vDV = -vDV;
%end

% Which side of the worm are we showing?
% view_angle = [0,90]; % XY
% view_angle = [0,0]; % XZ
% view_angle = [90,0]; % YZ
orient_vLR_text_halign = 'right';
orient_vLR_text_valign = 'bottom';
view_angle = [-25, 25]; % view the worm diagonally from above
if isHead
    view_angle = [25, 25]; % view the worm diagonally from above
end
if side_mode > 1
    
    % Which head neurons are we showing?
    iGanglia = [];
    if isHead
        iPharynx = contains({ganglion.name}, 'Pharyngeal');
        iDorsal = contains({ganglion.name}, 'Dorsal');
        iVentral = contains({ganglion.name}, 'Ventral Ganglion');
        iRVG = contains({ganglion.name}, 'Retrovesicular');
        if side_mode == 2
            iHead = contains({ganglion.name}, 'Anterior Ganglion (Left)') ...
                | contains({ganglion.name}, 'Lateral Ganglion (Left)');
        else
            iHead = contains({ganglion.name}, 'Anterior Ganglion (Right)') ...
                | contains({ganglion.name}, 'Lateral Ganglion (Right)');
        end
        iGanglia = iPharynx | iDorsal | iVentral | iRVG | iHead;
        
    % Which tail neurons are we showing?
    else
        iPAG = contains({ganglion.name}, 'Pre-Anal');
        iDRG = contains({ganglion.name}, 'Dorso-Rectal');
        if side_mode == 2
            iTail = contains({ganglion.name}, 'Tail (Left)');
        else
            iTail = contains({ganglion.name}, 'Tail (Right)');
        end
        iGanglia = iPAG | iDRG | iTail;
    end
    
    % Remove the neurons we're not showing.
    showNeurons = cat(1, ganglion(iGanglia).neurons);
    [~, ~, iShow] = intersect(showNeurons, neuron_names);
    neuron_names = neuron_names(iShow);
    neuron_sizes = neuron_sizes(iShow,:);
    neuron_RGBs = neuron_RGBs(iShow,:);
    neuron_alphas = neuron_alphas(iShow,:);
    gcamp = gcamp(iShow,:);
    pos = pos(iShow,:,:);
    
    % View the worm laterally.
    view_angle = [0, 0];
    
    % Set the left-right vector's text alignment.
    orient_vLR_text_halign = 'center';
    orient_vLR_text_valign = 'top';
end

% Translate the worm to have positive axes.
min_pos = min(min(pos, [], 3));
pos = pos - min_pos;

% Compute the the volume boundaries.
min_pos = min(min(pos, [], 3));
max_pos = max(max(pos, [], 3));
min_lim = min_pos - pos_lim_pad;
max_lim = max_pos + pos_lim_pad;

% Mirror the y-axis to match image conventions (vs. matrix conventions).
pos(:,2,:) = max_pos(2) - pos(:,2,:);

% Compute the orientation vector's for plotting.
orient_vec_size = 5;
orient_vec_pad = 3;
if isHead
    orient_vec_off = [max_pos(1) - min_pos(1) - orient_vec_pad; 0; ...
        max_pos(3) - orient_vec_size];
    orient_med_vAP = -vAP_med * orient_vec_size;
    orient_vAP = -vAP * orient_vec_size;
else
    orient_vec_off = [min_pos(1) + orient_vec_pad; 0; ...
        max_pos(3) - orient_vec_size];
    orient_med_vAP = vAP_med * orient_vec_size;
    orient_vAP = vAP * orient_vec_size;
end
orient_med_vLR = vLR_med * orient_vec_size;
orient_med_vDV = vDV_med * orient_vec_size;
orient_vLR = vLR * orient_vec_size;
orient_vDV = vDV * orient_vec_size;

% Organize the orientation text for plotting.
if isHead
    vAP_text = 'Anterior';
    orient_vAP_text = orient_med_vAP + orient_vec_off;
    orient_vAP_text_halign = 'right';
else
    vAP_text = 'Posterior';
    orient_vAP_text = orient_med_vAP + orient_vec_off;
    orient_vAP_text_halign = 'left';
end
vLR_text = 'Right';
vDV_text = 'Dorsal';
orient_vLR_text = orient_med_vLR + orient_vec_off;
orient_vDV_text = orient_med_vDV + orient_vec_off;

%% Organize the figure.
% Initialize the figure properties.
font_name = 'Arial';
font_size = 16;
neuron_font_size = 12;
if isHead
    neuron_font_size = 8;
end
figure_size = get(0, 'Screensize');
if ~isHead
    figure_size(3:4) = round(figure_size(3:4) .* 7/8);
    figure_size(1:2) = round(figure_size(1:2) + figure_size(3:4)./16);
end

% Are we making a movie?
movie = [];
if ~isempty(movie_name)
    movie = VideoWriter(movie_name, 'MPEG-4');
    movie.FrameRate = fps * (speed / compression);
    %movie.Quality = 10;
    open(movie);
end

% Which frames are we showing?
if isempty(frames)
    frames = 1:compression:size(pos,3);
    
% Check the bounding frames & compression.
else
    
    % Check the bounding frames.
    frames(1) = max(frames(1), 1);
    frames(end) = min(frames(end), size(pos,3));
    
    % Check the compression.
    if length(frames) > 1
        compression = round(mean(diff(frames)));
    end
end

% Create the figure.
if isempty(movie)
    figure('NumberTitle', 'off', 'Name', [exp_name ': ' plot_mode_label]);
    set(gcf, 'Position', figure_size);
else
    figure('Visible', 'off', ...
        'NumberTitle', 'off', 'Name', [exp_name ': ' plot_mode_label]);
    set(gcf, 'Position', [1 1 1920 1200]);
end
hold on;
axis equal;
xlim([min_lim(1), max_lim(1)]);
ylim([min_lim(2), max_lim(2)]);
zlim([min_lim(3), max_lim(3)]);
xlabel('X (microns)');
ylabel('Y (microns)');
zlabel('Z (microns)');
set(gca, 'FontName', font_name, 'FontSize', font_size);

%% Plot the experiment as a movie.
% Plot the neurons.
draw_speed = 1.0/fps;
for i=frames
    
    % Clear the figure.
    cla;
    
    % Is a stimulus being applied?
    stim_text = 'No Stimulus';
    stim_color = gray;
    for j = 1:length(stim)
        if i >= stim(j).on && i <= stim(j).off
            stim_text = stim_list(stim(j).num).name;
            stim_color = stim_list(stim(j).num).RGB;
        end
    end
    
    % How much time has elapsed?
    exp_time = i/fps;
    floor_exp_time = floor(exp_time);
    hours = floor(floor_exp_time/360);
    minutes = floor(floor_exp_time/60);
    seconds = mod(floor_exp_time,60);
    milliseconds = round((exp_time - floor_exp_time) * 1000);
    
    % Show the stimulus and time.
    title_text = sprintf('%s Neurons (Frame=%04d Time=%02d:%02d:%02d:%03d) - %s', ...
        body_text, i, hours, minutes, seconds, milliseconds, upper(stim_text));
    title(title_text, 'FontWeight', 'bold');
    set(gcf, 'Color', stim_color);

    % Plot the median orientation.
    quiver3(orient_vec_off(1), orient_vec_off(2), orient_vec_off(3), ...
        orient_med_vAP(1), orient_med_vAP(2), orient_med_vAP(3), ...
        'LineStyle', ':', 'LineWidth', 3, 'Color', 'r', 'MaxHeadSize', 1);
    quiver3(orient_vec_off(1), orient_vec_off(2), orient_vec_off(3), ...
        orient_med_vLR(1), orient_med_vLR(2), orient_med_vLR(3), ...
        'LineStyle', ':', 'LineWidth', 3, 'Color', 'g', 'MaxHeadSize', 1);
    quiver3(orient_vec_off(1), orient_vec_off(2), orient_vec_off(3), ...
        orient_med_vDV(1), orient_med_vDV(2), orient_med_vDV(3), ...
        'LineStyle', ':', 'LineWidth', 3, 'Color', 'b', 'MaxHeadSize', 1);
    
    % Plot the orientation.
    quiver3(orient_vec_off(1), orient_vec_off(2), orient_vec_off(3), ...
        orient_vAP(1,i), orient_vAP(2,i), orient_vAP(3,i), ...
        'LineWidth', 5, 'Color', 'r', 'MaxHeadSize', 1);
    quiver3(orient_vec_off(1), orient_vec_off(2), orient_vec_off(3), ...
        orient_vLR(1,i), orient_vLR(2,i), orient_vLR(3,i), ...
        'LineWidth', 5, 'Color', 'g', 'MaxHeadSize', 1);
    quiver3(orient_vec_off(1), orient_vec_off(2), orient_vec_off(3), ...
        orient_vDV(1,i), orient_vDV(2,i), orient_vDV(3,i), ...
        'LineWidth', 5, 'Color', 'b', 'MaxHeadSize', 1);
    
    % Label the orientation.
    text(orient_vAP_text(1), orient_vAP_text(2), orient_vAP_text(3), ...
        vAP_text, 'HorizontalAlignment', orient_vAP_text_halign, ...
         'Color', 'r', 'FontName', font_name, 'FontSize', font_size);
    text(orient_vLR_text(1), orient_vLR_text(2), orient_vLR_text(3), ...
        vLR_text, 'HorizontalAlignment', orient_vLR_text_halign, ...
        'VerticalAlignment', orient_vLR_text_valign, ...
        'Color', 'g', 'FontName', font_name, 'FontSize', font_size);
    text(orient_vDV_text(1), orient_vDV_text(2), orient_vDV_text(3), ...
        vDV_text, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Color', 'b', 'FontName', font_name, 'FontSize', font_size);
    
    % Plot the neurons.
    for j = 1:size(pos,1)
        if  ~isnan(gcamp(j,i)) && ~any(isnan(pos(j,:,i)))
            patch('Faces', neuron_patch.faces, 'Vertices', ...
                neuron_patch.vertices * neuron_sizes(j,i) + pos(j,:,i), ...
                'FaceColor', neuron_RGBs(j,:), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', neuron_alphas(j,i));
            text(pos(j,1,i), pos(j,2,i), pos(j,3,i), neuron_names{j}, ...
                'FontWeight', 'bold', ...
                'FontSize', neuron_font_size, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end
    
    % Shade the figure.
    view(view_angle);
    camlight;
    lighting phong; %gouraud;
    material shiny; %dull;
    
    % Draw everything.
    if isempty(movie)
        drawnow;
        pause(draw_speed);
    else
        writeVideo(movie, getframe(gcf));
        if mod(i,100) == 0
            disp(['Writing frame #: ' num2str(i) '/' num2str(length(frames))]);
        end
    end
end

% Close the movie.
if ~isempty(movie)
    close(movie);
end
end