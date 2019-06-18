%% Initialize the data.

% Organize the trace files.
files = natsort(split(ls('*_traces.mat')));
files = files(~cellfun('isempty', files));

% Are we smoothing?
is_smooth = true;

% Load the data.
if is_smooth
    load('exp_data_smooth.mat');
else
    load('exp_data.mat');
end

% Set the passband.
passband = false;
if is_smooth
    passband = [2, 240, 1];
end

% Translate experiments to indices.
exps_i = nan(66,1);
exps_i(1) = 1;
exps_i(22) = 2;
exps_i(23) = 3;
exps_i(55) = 4;
exps_i(56) = 5;
exps_i(65) = 6;
exps_i(66) = 7;


%% Plot locomotion neurons.
% Head #55:
% OLQVL/~R/~D,
% ~RMEL/R/~D/V
% RID
% AVA
% AVB
% AVD~L/R (remove light ~5s)
% AVE
% RIM
% ~RIR
% ~RIV (remove light ~5s)
% RMDL/R
% SMDV
% SAAV~R
% SAAD~L
% SIADR/VL
% SMDD
% DA1
% DB1
% DB2
% RIF
% SABV/D
% VB1
% VB2

% Choose the neurons/experiment.
names = {
    'OLQ'
    'RME'
    'RID'
    'AVA'
    'AVB'
    'AVD'
    'AVE'
    'RIM'
    'RIR'
    'RIV'
    'RMD'
    'SMD'
    'SAA'
    'SIA'
    'DA1'
    'DB1'
    'DB2'
    'RIF'
    'SAB'
    'VB1'
    'VB2'
    };
exp_name = '55_head_run301_traces.mat';

% Plot the data.
plotExpTraces(exp_name, false, true, 3, passband, names, 10);
box on
ticks = 20:30:230;
tick_labels = arrayfun(@num2str, ticks + 10, 'UniformOutput', false);
xticks(ticks - 0.1);
xticklabels(tick_labels);
%title('Left/Right Locomotion');

% Choose a smaller subset of the neurons.
names = {
    'AVAR'
    'AVBR'
    'AVDR'
    'AVER'
    'DA1' %~
    'DB1' %~
    'DB2' %~
    'OLQVL'
    'RID'
    'RIFL' %~
    'RIMR'
    %'RIR'
    'RIVR'
    'RMDR'
    'RMER'
    'RMEV'
    'SAAVR' %~
    'SABD'
    'SABVR'
    'SIADR'
    'SMDDL'
    'SMDVL'
    'VB1'
    'VB2'
    };

% Plot the data.
plotExpTraces(exp_name, false, true, 3, passband, names, 10);
box on
ticks = 20:30:230;
tick_labels = arrayfun(@num2str, ticks + 10, 'UniformOutput', false);
xticks(ticks - 0.1);
xticklabels(tick_labels);
%title('Locomotion');

%% Plot NaCl Exemplary Responses.
% ASEL: 55/56L
% X ASER: 23/55/66R
% X ASE: *56R
% ADF: 66R
% X ~AIN: 56R
% ASG: 55R
% ASI: 66R
% ASK: 56L
% AUA: 55L
% AVH: 55/66L
% AVJ: 55R
% AWB: 65L/66R
% ~IL2D/V: 22R/66R
% RIB: 55R

% Choose the neurons/experiments.
names = {
    'ASEL'
    'ADFR'
    'ASGR'
    'ASIR'
    'ASKL'
    'AUAL'
    'AVHL'
    'AVJR'
    'AWBR'
    'IL2DR'
    'IL2VR'
    'RIBR'
    };
exp_nums = [
    56
    66
    55
    66
    56
    55
    55
    55
    65
    22
    66
    55
    ];
stims = ones(size(exp_nums)) * 3;

% Plot the data.
plotAnimalResponses(neurons, names, exps_i(exp_nums), stims, fps, true);
pos = get(gcf, 'Position');
pos(3) = pos(3) / 2;
set(gcf, 'Position', pos);
title('NaCl Exemplary Responses');
box on;


%% Plot Published AWC Repsonses.
% AWC(but): 65R
% AWC(pent): 65L

% Choose the neurons/experiments.
names = {
%    'AWCL'
    'AWCR'
    'AWCL'
%    'AWCR'
    };
exp_nums = [
%    65
    65
    65
%    65
    ];
stims = [
%    1
    1
    2
%    2
    ];

% Plot the data.
plotAnimalResponses(neurons, names, exps_i(exp_nums), stims, fps, true);
pos = get(gcf, 'Position');
pos(3) = pos(3) / 3;
pos(4) = pos(4) / 4;
set(gcf, 'Position', pos);
title('Published AWC Responses');
box on;


%% Plot Published AWA Repsonses.
% AWA(but): 56R(smooth)
% AWA(but): 65R(spiky)
% AWA(pent): 56/65/66R

% Choose the neurons/experiments.
names = {
    'AWAR'
    'AWAR'
    'AWAR'
    };
exp_nums = [
    65
    56
    56
    ];
stims = [
    1
    1
    2
    ];

% Plot the data.
plotAnimalResponses(neurons, names, exps_i(exp_nums), stims, fps, true);
pos = get(gcf, 'Position');
pos(3) = pos(3) / 3;
pos(4) = pos(4) / 3;
set(gcf, 'Position', pos);
title('Published AWA Responses');
box on;


%% Plot non-significant interesting neurons.
% ADL(NaCl): 55L
% ADL(pent): 55L
% M2L(NaCl): 66L
% X ASH(NaCl): 55L
% X ASH(pent): 55L

% Choose the neurons/experiments.
names = {
    'ADLL'
    'ADLL'
    'M2L'
    };
exp_nums = [
    55
    55
    66
    ];
stims = [
    3
    2
    3
    ];

% Plot the data.
plotAnimalResponses(neurons, names, exps_i(exp_nums), stims, fps, true);
pos = get(gcf, 'Position');
pos(3) = pos(3) / 2;
pos(4) = pos(4) / 3;
set(gcf, 'Position', pos);
title('Non-significant Interesting Responses');
box on;


%% Plot spiking neurons.
% I5 (vs. earlier AWA): 66(pent to end) or 65(3 spikes)

% Choose the neurons/experiment.
names = {
    'AWAL'
    'I5'
    };
exp_name = '66_head_run101_traces.mat';

% Plot the data.
plotExpTraces(exp_name, false, true, 3, passband, names, 80);
pos = get(gcf, 'Position');
pos(4) = pos(4) / 3;
set(gcf, 'Position', pos);
title('Spiking Neurons');
box on;
ticks = 20:30:160;
tick_labels = arrayfun(@num2str, ticks + 70, 'UniformOutput', false);
xticks(ticks - 1/3);
xticklabels(tick_labels);

% Plot the raw data.
plotExpTraces(exp_name, false, true, 3, false, names, 80);
pos = get(gcf, 'Position');
pos(4) = pos(4) / 3;
set(gcf, 'Position', pos);
title('Spiking Neurons');
box on;
ticks = 20:30:160;
tick_labels = arrayfun(@num2str, ticks + 70, 'UniformOutput', false);
xticks(ticks - 1/3);
xticklabels(tick_labels);

%% Plot pharyngeal responses to light.
% I6: unresponsive
% I2: light 56
% M1: light 55, 56, 65
% MI: light 22, 55, 56

% Choose the neurons/experiment.
names = {
    %'ADF'
    %'ADLL'
    %'ASIL'
    'NSM'
    'I2R'
    'M1'
    'MI'
    };
exp_name = '56_head_run401_traces.mat';

% Plot the data.
plotExpTraces(exp_name, false, true, 3, passband, names);
pos = get(gcf, 'Position');
pos(4) = pos(4) / 2;
set(gcf, 'Position', pos);
%title('Light Responses');
box on;
xticks(0:30:240);
