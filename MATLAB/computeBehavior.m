function behavior = computeBehavior(neurons, traces, positions, fps, ...
    varargin)
%COMPUTEBEHAVIOR Compute a worm's behavior.
%
%   BEHAVIORS = COMPUTEBEHAVIOR(NEURONS, TRACES, POSITIONS, FPS)
%   BEHAVIORS = COMPUTEBEHAVIOR(NEURONS, TRACES, POSITIONS, FPS,
%                               ISCLEANTRACES)
%
%   Inputs:
%   neurons       = a cell array of neuron names
%   traces        = the matrix of neural activity
%   positions     = the matrix of neural positions
%   fps           = frames/second
%   isCleanTraces = should we clean the traces?
%                   default = true
%
%   Output:
%   behavior = a struct array of behaviors
%       name   = the name of the behavior
%       ylabel = the y-axis label for the behavior
%       data   = the data for the behavior

% Should we clean the traces?
isCleanTraces = true;
if ~isempty(varargin)
    isCleanTraces = varargin{1};
end

% Initialize the behavior & neurons.
behavior = [];
if iscategorical(neurons)
    neurons = cellstr(neurons);
end
neurons_LR = neurons;
neurons = stripNeuronLR(neurons, false);

% Clean up the neural traces.
if isCleanTraces
    traces = cleanTraces(traces);
    [iNeuron, iFrame] = find(isnan(traces));
    for i = 1:length(iNeuron)
        positions(iNeuron(i),:,iFrame(i)) = nan;
    end
end
i = 1;

% Compute the global activity.
behavior(i).name = 'Global Brain Activity';
behavior(i).ylabel = 'Neural Activity (scaled fluoresence)';
behavior(i).data = rescale(nanmean(traces, 1));
behavior(i).center = nanmedian(behavior(i).data);
i = i + 1;

% Compute the global x-axis location.
behavior(i).name = 'Global X-Axis';
behavior(i).ylabel = 'X (microns)';
behavior(i).data = nanmean(squeeze(positions(:,1,:)), 1);
behavior(i).center = nanmedian(behavior(i).data);
i = i + 1;

% Compute the global y-axis location.
behavior(i).name = 'Global Y-Axis';
behavior(i).ylabel = 'Y (microns)';
behavior(i).data = nanmean(squeeze(positions(:,2,:)), 1);
behavior(i).center = nanmedian(behavior(i).data);
i = i + 1;

% Compute the global z-axis location.
behavior(i).name = 'Global Z-Axis';
behavior(i).ylabel = 'Z (microns)';
behavior(i).data = nanmean(squeeze(positions(:,3,:)), 1);
behavior(i).center = nanmedian(behavior(i).data);
i = i + 1;

% Compute the global motion.
behavior(i).name = 'Global Motion';
behavior(i).ylabel = '\deltaNeuronal Positions/Seconds (microns)';
round_fps = round(fps);
mean_pos = squeeze(nanmean(positions, 1));
diff_pos = diff(mean_pos, round_fps, 2);
behavior(i).data = cat(2, nan(1, round_fps), sqrt(nansum(diff_pos.^2, 1)));
behavior(i).center = nanmedian(behavior(i).data);
i = i + 1;

% Compute the anterior pharyngeal expansion & contraction.
%ant_pharynx = {'I2L', 'I2R', 'I3','M3L', 'M3R'};
ant_pharynx = {'M3L', 'M3R'};
[~, ~, iAP] = intersect(ant_pharynx, neurons_LR);
if length(iAP) > 1
    behavior(i).name = 'Anterior Pharyngeal Expansion (+) & Contraction (-)';
    behavior(i).ylabel = 'Size (microns)';
    data = computeDistance(positions(iAP,:,:));
    behavior(i).data = data - nanmedian(data);
    behavior(i).center = 0;
    i = i + 1;
end

% Compute the posterior pharyngeal expansion & contraction.
%post_pharynx = {'I6','M1','M5'};
post_pharynx = {'M1','M5'};
[~, ~, iPP] = intersect(post_pharynx, neurons);
if length(iPP) > 1
    behavior(i).name = 'Posterior Pharyngeal Expansion (+) & Contraction (-)';
    behavior(i).ylabel = 'Size (microns)';
    data = computeDistance(positions(iPP,:,:));
    behavior(i).data = data - nanmedian(data);
    behavior(i).center = nanmedian(data);
    i = i + 1;
end

% Compute the head's posture.
ant_head_L = {'OLQDL','OLQVL'};
ant_head_R = {'OLQDR','OLQVR'};
post_head_L = {'AVAL','AVEL'};
post_head_R = {'AVAR','AVER'};
[~, ~, iAHL] = intersect(ant_head_L, neurons_LR);
[~, ~, iAHR] = intersect(ant_head_R, neurons_LR);
[~, ~, iPHL] = intersect(post_head_L, neurons_LR);
[~, ~, iPHR] = intersect(post_head_R, neurons_LR);
is_head_pos = false;
if length(iAHL) == length(ant_head_L) && ...
        length(iAHR) == length(ant_head_R) && ...
        length(iPHL) == length(post_head_L) && ...
        length(iPHR) == length(post_head_R)
    is_head_pos = true;
    
    % Compute the head neurons' positons.
    pos_AHL = squeeze(mean(positions(iAHL,:,:), 1));
    pos_AHR = squeeze(mean(positions(iAHR,:,:), 1));
    pos_PHL = squeeze(mean(positions(iPHL,:,:), 1));
    pos_PHR = squeeze(mean(positions(iPHR,:,:), 1));
        
    % Compute the vectors for the worm's orientation.
    vAPL = pos_PHL - pos_AHL; % AP vectors on left
    vAPR = pos_PHR - pos_AHR; % AP vectors on right
    vAP = (vAPL + vAPR) / 2; % AP vectors
    vAP = vAP ./ sqrt(nansum(vAP.^2,1)); % AP unit vectors

    % Compute the vectors for the worm's L-R orientation.
    vLR = pos_PHR - pos_PHL; % LR vector
       
    % Orient LR orthogonal to AP by removing its AP projection.
    vLR = vLR - squeeze(dot(vLR, vAP)) .* vAP; % LR vectors
    vLR = vLR ./ sqrt(sum(vLR.^2,1)); % LR unit vectors
    
    % Orient DV orthogonal to AP and LR.
    vDV = cross(vLR, vAP); % DV vectors
    
    % Compute the median orientation unit vectors.
    %vAP_med = median(vAP,2);
    %vLR_med = median(vLR,2);
    %vLR_med = vLR_med - squeeze(dot(vLR_med, vAP_med)) .* vAP_med; % LR vectors
    %vDV_med = cross(vLR_med, vAP_med);
    %vAP_med = vAP_med ./ sqrt(sum(vAP_med.^2,1));
    %vLR_med = vLR_med ./ sqrt(sum(vLR_med.^2,1));
    %vDV_med = vDV_med ./ sqrt(sum(vDV_med.^2,1));
    vAP_med = vAP(:,nanMedianIndex(vAP,2));
    vLR_med = vLR(:,nanMedianIndex(vLR,2));
    
    % Orient LR orthogonal to AP by removing its AP projection.
    vLR_med = vLR_med - dot(vLR_med, vAP_med) .* vAP_med;
    vLR_med = vLR_med ./ sqrt(sum(vLR_med.^2, 1));
    
    % Orient DV orthogonal to AP and LR.
    vDV_med = cross(vLR_med, vAP_med);
    
    % Create an array of the median vectors for later computations.
    vAP_meds = repmat(vAP_med, 1, size(vAP,2));
    vDV_meds = repmat(vDV_med, 1, size(vAP,2));
    vLR_meds = repmat(vLR_med, 1, size(vAP,2));
end

% Compute head extension & retraction.
% Head extension & retraction is defined as the mean of the distance,
% between the anterior & posterior head neurons, on both sides of the worm,
% centered at the median of this measure.
if is_head_pos
    
    % Compute the distance between the anterior & posterior of the head.
    dist_HL = sqrt(sum((pos_AHL - pos_PHL).^2, 1));
    dist_HR = sqrt(sum((pos_AHR - pos_PHR).^2, 1));    
    dist_H = (dist_HL + dist_HR) ./ 2;
    med_dist_H = nanmedian(dist_H);
    data = dist_H - med_dist_H;
    center = med_dist_H;
    
    % Compute head extension & retraction.
    behavior(i).name = 'Head Extension (+) & Retraction (-)';
    behavior(i).ylabel = 'A-P Distance from Neutral (microns)';
    behavior(i).data = data;
    behavior(i).center = center;
    i = i + 1;
end

% Compute head twist, its rotation in the anterior-posterior plane (x axis).
% Head twist is defined as the angle (in radians), between the worm's
% DV vector (z axis) and the median of this vector, projected in the
% LR-DV plane (z in y-z axis).
if is_head_pos
    
    % Project the DV vector into the AP plane.
    pDV = vDV - dot(vDV, vAP_meds) .* vAP_meds;
    
    % Compute the DV vector angle in the LR-DV plane (z in y-z axis).
    cross_pDV = cross(pDV, vDV_meds);
    mag_cross_pDV = sqrt(sum(cross_pDV.^2, 1));
    angle_sign = sign(dot(cross_pDV, vAP_meds));
    angles = atan2(mag_cross_pDV, dot(pDV, vDV_meds));
    data = angles .* angle_sign;
    center = nanmedian(data);
    
    % Compute the head twist.
    behavior(i).name = 'Head Twist: Up Left (-) & Up Right (+)';
    behavior(i).ylabel = 'A-P Angle from Neutral (radians)';
    behavior(i).data = data;
    behavior(i).center = center;
    i = i + 1;
end

% Compute dorsal-ventral head motion.
% Dorsal-ventral head motion is defined as the angle (in radians), between
% the worm's AP vector (x axis) and  the median of this vector, projected
% in the AP-DV plane (x in x-z axis).
if is_head_pos
    
    % Project the AP vector into the LR plane.
    pAP = vAP  - dot(vAP, vLR_meds) .* vLR_meds;
    
    % Compute the AP vector angle in the AP-DV plane (x in x-z axis).
    cross_pAP = cross(pAP, vAP_meds);
    mag_cross_pAP = sqrt(sum(cross_pAP.^2, 1));
    angle_sign = sign(dot(cross_pAP, vLR_meds));
    angles = atan2(mag_cross_pAP, dot(pAP, vAP_meds));
    data = angles .* angle_sign;
    center = nanmedian(data);
    
    % Compute the tail twist.
    behavior(i).name = 'Head Dorsal (+) & Ventral (-) Motion';
    behavior(i).ylabel = 'D-V Angle from Neutral (radians)';
    behavior(i).data = data;
    behavior(i).center = center;
    i = i + 1;
end

% Compute left-right head motion.
% Left-right head motion is defined as the angle (in radians), between
% the worm's LR vector (y axis) and the median of this vector, projected
% in the AP-LR plane (y in x-y axis).
if is_head_pos
    
    % Project the LR vector into the DV plane.
    pLR = vLR  - dot(vLR, vDV_meds) .* vDV_meds;
    
    % Compute the LR vector angle in the AP-LR plane (y in x-y axis).
    cross_pLR = cross(pLR, vLR_meds);
    mag_cross_pLR = sqrt(sum(cross_pLR.^2, 1));
    angle_sign = -sign(dot(cross_pLR, vDV_meds));
    angles = atan2(mag_cross_pLR, dot(pLR, vLR_meds));
    data = angles .* angle_sign;
    center = nanmedian(data);
    
    % Compute the tail twist.
    behavior(i).name = 'Head Left (-) & Right (+) Motion';
    behavior(i).ylabel = 'L-R Angle from Neutral (radians)';
    behavior(i).data = data;
    behavior(i).center = center;
end

% Compute the tail's posture.
ant_tail_L = {'LUAL','PVCL'};
ant_tail_R = {'LUAR','PVCR'};
post_tail_L = {'PVNL','PLML'};
post_tail_R = {'PVNR','PLMR'};
[~, ~, iATL] = intersect(ant_tail_L, neurons_LR);
[~, ~, iATR] = intersect(ant_tail_R, neurons_LR);
[~, ~, iPTL] = intersect(post_tail_L, neurons_LR);
[~, ~, iPTR] = intersect(post_tail_R, neurons_LR);
is_tail_pos = false;
if length(iATL) == length(ant_tail_L) && ...
        length(iATR) == length(ant_tail_R) && ...
        length(iPTL) == length(post_tail_L) && ...
        length(iPTR) == length(post_tail_R)
    is_tail_pos = true;
    
    % Compute the tail neurons' positons.
    pos_ATL = squeeze(mean(positions(iATL,:,:), 1));
    pos_ATR = squeeze(mean(positions(iATR,:,:), 1));
    pos_PTL = squeeze(mean(positions(iPTL,:,:), 1));
    pos_PTR = squeeze(mean(positions(iPTR,:,:), 1));
    
    % Compute the vectors for the worm's A-P orientation.
    vAPL = pos_PTL - pos_ATL; % AP vectors on left
    vAPR = pos_PTR - pos_ATR; % AP vectors on right
    vAP = (vAPL + vAPR) / 2; % AP vectors
    vAP = vAP ./ sqrt(sum(vAP.^2,1)); % AP unit vectors
    
    % Compute the vectors for the worm's L-R orientation.
    vLR = pos_ATR - pos_ATL; % LR vector
       
    % Orient LR orthogonal to AP by removing its AP projection.
    vLR = vLR - squeeze(dot(vLR, vAP)) .* vAP; % LR vectors
    vLR = vLR ./ sqrt(sum(vLR.^2,1)); % LR unit vectors
    
    % Orient DV orthogonal to AP and LR.
    vDV = cross(vLR, vAP); % DV vectors
    
    % Compute the median orientation unit vectors.
    %vAP_med = median(vAP,2);
    %vLR_med = median(vLR,2);
    %vLR_med = vLR_med - squeeze(dot(vLR_med, vAP_med)) .* vAP_med; % LR vectors
    %vDV_med = cross(vLR_med, vAP_med);
    %vAP_med = vAP_med ./ sqrt(sum(vAP_med.^2,1));
    %vLR_med = vLR_med ./ sqrt(sum(vLR_med.^2,1));
    %vDV_med = vDV_med ./ sqrt(sum(vDV_med.^2,1));
    vAP_med = vAP(:,nanMedianIndex(vAP,2));
    vLR_med = vLR(:,nanMedianIndex(vLR,2));
    
    % Orient LR orthogonal to AP by removing its AP projection.
    vLR_med = vLR_med - dot(vLR_med, vAP_med) .* vAP_med;
    vLR_med = vLR_med ./ sqrt(sum(vLR_med.^2, 1));
    
    % Orient DV orthogonal to AP and LR.
    vDV_med = cross(vLR_med, vAP_med);
    
    % Create an array of the median vectors for later computations.
    vAP_meds = repmat(vAP_med, 1, size(vAP,2));
    vDV_meds = repmat(vDV_med, 1, size(vAP,2));
    vLR_meds = repmat(vLR_med, 1, size(vAP,2));
end

% Compute tail extension & retraction.
% Tail extension & retraction is defined as the mean of the distance,
% between the anterior & posterior tail neurons, on both sides of the worm,
% centered at the median of this measure.
if is_tail_pos
    
    % Compute the distance between the anterior & posterior of the tail.
    dist_TL = sqrt(sum((pos_ATL - pos_PTL).^2, 1));
    dist_TR = sqrt(sum((pos_ATR - pos_PTR).^2, 1));    
    dist_T = (dist_TL + dist_TR) ./ 2;
    med_dist_T = nanmedian(dist_T);
    data = dist_T - med_dist_T;
    center = med_dist_T;
    
    % Compute tail extension & retraction.
    behavior(i).name = 'Tail Extension (+) & Retraction (-)';
    behavior(i).ylabel = 'A-P Distance from Neutral (microns)';
    behavior(i).data = data;
    behavior(i).center = center;
    i= i + 1;
end

% Compute tail twist, its rotation in the anterior-posterior plane (x axis).
% Tail twist is defined as the angle (in radians), between the worm's
% DV vector (z axis) and the median of this vector, projected in the
% LR-DV plane (z in y-z axis).
if is_tail_pos
    
    % Project the DV vector into the AP plane.
    pDV = vDV - dot(vDV, vAP_meds) .* vAP_meds;
    
    % Compute the DV vector angle in the LR-DV plane (z in y-z axis).
    cross_pDV = cross(pDV, vDV_meds);
    mag_cross_pDV = sqrt(sum(cross_pDV.^2, 1));
    angle_sign = sign(dot(cross_pDV, vAP_meds));
    angles = atan2(mag_cross_pDV, dot(pDV, vDV_meds));
    data = angles .* angle_sign;
    center = nanmedian(data);
    
    % Compute the tail twist.
    behavior(i).name = 'Tail Twist: Up Left (-) & Up Right (+)';
    behavior(i).ylabel = 'A-P Angle from Neutral (radians)';
    behavior(i).data = data;
    behavior(i).center = center;
    i = i + 1;
end

% Compute dorsal-ventral tail motion.
% Dorsal-ventral tail motion is defined as the angle (in radians), between
% the worm's AP vector (x axis) and  the median of this vector, projected
% in the AP-DV plane (x in x-z axis).
if is_tail_pos
    
    % Project the AP vector into the LR plane.
    pAP = vAP  - dot(vAP, vLR_meds) .* vLR_meds;
    
    % Compute the AP vector angle in the AP-DV plane (x in x-z axis).
    cross_pAP = cross(pAP, vAP_meds);
    mag_cross_pAP = sqrt(sum(cross_pAP.^2, 1));
    angle_sign = -sign(dot(cross_pAP, vLR_meds));
    angles = atan2(mag_cross_pAP, dot(pAP, vAP_meds));
    data = angles .* angle_sign;
    center = nanmedian(data);
    
    % Compute the tail twist.
    behavior(i).name = 'Tail Dorsal (+) & Ventral (-) Motion';
    behavior(i).ylabel = 'D-V Angle from Neutral (radians)';
    behavior(i).data = data;
    behavior(i).center = center;
    i = i + 1;
end

% Compute left-right tail motion.
% Left-right tail motion is defined as the angle (in radians), between
% the worm's LR vector (y axis) and the median of this vector, projected
% in the AP-LR plane (y in x-y axis).
if is_tail_pos
    
    % Project the LR vector into the DV plane.
    pLR = vLR  - dot(vLR, vDV_meds) .* vDV_meds;
    
    % Compute the LR vector angle in the AP-LR plane (y in x-y axis).
    cross_pLR = cross(pLR, vLR_meds);
    mag_cross_pLR = sqrt(sum(cross_pLR.^2, 1));
    angle_sign = sign(dot(cross_pLR, vDV_meds));
    angles = atan2(mag_cross_pLR, dot(pLR, vLR_meds));
    data = angles .* angle_sign;
    center = nanmedian(data);
    
    % Compute the tail twist.
    behavior(i).name = 'Tail Left (-) & Right (+) Motion';
    behavior(i).ylabel = 'L-R Angle from Neutral (radians)';
    behavior(i).data = data;
    behavior(i).center = center;
end
end

%% Compute the mean distance between all neurons.
function dist = computeDistance(positions)
dist = nan(1, size(positions,3));
switch size(positions,1)
    
    % Not enough neurons.
    case 1
        % do nothing

    % Compute the distance between the two neurons.
    case 2
        pos_diff = positions(1,:,:) - positions(2,:,:);
        dist(1,:) = sqrt(squeeze(sum(pos_diff .^ 2, 2)));
        
    % Compute the mean distance between all neurons.
    otherwise
        dist = nan(size(positions,1) - 1, size(positions,3));
        for i = 1:size(dist,1)
            pos_diff = positions(i,:,:) - positions((i+1):end,:,:);
            pos_dist = sqrt(sum(pos_diff .^ 2, 2));
            dist(i,:) = nanmean(pos_dist, 1);
        end
        dist = nanmean(dist, 1);
end
end

%% Rescale data to [0,1].
function data = rescale(data)
data = data - nanmin(data);
data = data ./ nanmax(data);
end
