% Load the stats.
load('exp_stats.mat');


%% Print the neuron stats.

% Create the neuron stats table.
stats_table = table;
stats_table.Neuron = categorical(stats.neuron);

% Determine the significant neurons.
but_sig = stats.on.qT.on(:,1) < 0.1;
pent_sig = stats.on.qT.on(:,2) < 0.1;
nacl_sig = stats.on.qT.on(:,3) < 0.1;
light_sig = stats.on.qT.on(:,4) < 0.1;

% Butanone stats table.
stats_table.Butanone_Samples = stats.N(:,1);
stats_table.Butanone_Mean = stats.on.mean(:,1);
stats_table.Butanone_Min_CI = stats.on.ci(:,1,1);
stats_table.Butanone_Max_CI = stats.on.ci(:,1,2);
stats_table.Butanone_p = stats.on.pT(:,1);
stats_table.Butanone_q = stats.on.qT.on(:,1);
stats_table.Butanone_Sig = but_sig;

% Pentanedione stats table.
stats_table.Pentanedione_Samples = stats.N(:,2);
stats_table.Pentanedione_Mean = stats.on.mean(:,2);
stats_table.Pentanedione_Min_CI = stats.on.ci(:,2,1);
stats_table.Pentanedione_Max_CI = stats.on.ci(:,2,2);
stats_table.Pentanedione_p = stats.on.pT(:,2);
stats_table.Pentanedione_q = stats.on.qT.on(:,2);
stats_table.Pentanedione_Sig = pent_sig;

% NaCl stats table.
stats_table.NaCl_Samples = stats.N(:,3);
stats_table.NaCl_Mean = stats.on.mean(:,3);
stats_table.NaCl_Min_CI = stats.on.ci(:,3,1);
stats_table.NaCl_Max_CI = stats.on.ci(:,3,2);
stats_table.NaCl_p = stats.on.pT(:,3);
stats_table.NaCl_q = stats.on.qT.on(:,3);
stats_table.NaCl_Sig = nacl_sig;

% Light stats table.
stats_table.Light_Samples = stats.N(:,4);
stats_table.Light_Mean = stats.on.mean(:,4);
stats_table.Light_Min_CI = stats.on.ci(:,4,1);
stats_table.Light_Max_CI = stats.on.ci(:,4,2);
stats_table.Light_p = stats.on.pT(:,4);
stats_table.Light_q = stats.on.qT.on(:,4);
stats_table.Light_Sig = light_sig;

% Print the stats tables.
writetable(stats_table, 'stim_stats.csv');


%% Print the significant neuron's sign.

% Create the neuron signs table.
signs_table = table;
signs_table.Neuron = categorical(stats.neuron);

% Initialize the signs.
no_sign = categorical(" ");
neg_sign = categorical("-");
pos_sign = categorical("+");

% Butanone signs table.
signs_table.Butanone(:) = no_sign;
but_neg = but_sig & stats.on.ci(:,1) < 0;
but_pos = but_sig & stats.on.ci(:,1) > 0;
signs_table.Butanone(but_neg) = neg_sign;
signs_table.Butanone(but_pos) = pos_sign;

% Pentanedione signs table.
signs_table.Pentanedione(:) = no_sign;
pent_neg = pent_sig & stats.on.ci(:,2) < 0;
pent_pos = pent_sig & stats.on.ci(:,2) > 0;
signs_table.Pentanedione(pent_neg) = neg_sign;
signs_table.Pentanedione(pent_pos) = pos_sign;

% NaCl signs table.
signs_table.NaCl(:) = no_sign;
nacl_neg = nacl_sig & stats.on.ci(:,3) < 0;
nacl_pos = nacl_sig & stats.on.ci(:,3) > 0;
signs_table.NaCl(nacl_neg) = neg_sign;
signs_table.NaCl(nacl_pos) = pos_sign;

% Light signs table.
signs_table.Light(:) = no_sign;
light_neg = light_sig & stats.on.ci(:,4) < 0;
light_pos = light_sig & stats.on.ci(:,4) > 0;
signs_table.Light(light_neg) = neg_sign;
signs_table.Light(light_pos) = pos_sign;

% Remove insignificant neurons from the table.
signs_table = signs_table(but_sig | pent_sig | nacl_sig,:);

% Print the signs tables.
writetable(signs_table, 'stim_signs.csv');
