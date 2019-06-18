%% Show the significant neurons & their confidence intervals.
load('exp_stats.mat');

% q-significance threshold.
q_thresh = 0.1;
disp(['*** Q = ' num2str(q_thresh)]);

% All neurons.
neuron_on = zeros(size(stats.neuron));
neuron_off = zeros(size(stats.neuron));

% Butanone ON.
stim = 1;
neuron_i = find(stats.on.qT.on(:,stim) <= q_thresh);
neuron_on(neuron_i) = neuron_on(neuron_i) + 1;
disp(['*** 2-butanone ON = ' num2str(length(neuron_i))]);
for i = 1:length(neuron_i)
    j = neuron_i(i);
    str = sprintf('%s (q=%.3f, ci=%d [%.2f,%.2f], N=%d)', ...
        stats.neuron{j}, stats.on.qT.on(j,stim), ...
        stats.on.ci_sign(j,stim), ...
        stats.on.ci(j,stim,1), stats.on.ci(j,stim,2), stats.N(j,stim));
    disp(str);
end
disp(' ');

% Butanone OFF.
stim = 1;
neuron_i = find(stats.off.qT.off(:,stim) <= q_thresh);
neuron_off(neuron_i) = neuron_off(neuron_i) + 1;
disp(['*** 2-butanone OFF = ' num2str(length(neuron_i))]);
for i = 1:length(neuron_i)
    j = neuron_i(i);
    str = sprintf('%s (q=%.3f, ci=%d [%.2f,%.2f], N=%d)', ...
        stats.neuron{j}, stats.off.qT.off(j,stim), ...
        stats.off.ci_sign(j,stim), ...
        stats.off.ci(j,stim,1), stats.off.ci(j,stim,2), stats.N(j,stim));
    disp(str);
end
disp(' ');

% Pentanedione ON.
stim = 2;
neuron_i = find(stats.on.qT.on(:,stim) <= q_thresh);
neuron_on(neuron_i) = neuron_on(neuron_i) + 1;
disp(['*** 2,3-pentanedione ON = ' num2str(length(neuron_i))]);
for i = 1:length(neuron_i)
    j = neuron_i(i);
    str = sprintf('%s (q=%.3f, ci=%d [%.2f,%.2f], N=%d)', ...
        stats.neuron{j}, stats.on.qT.on(j,stim), ...
        stats.on.ci_sign(j,stim), ...
        stats.on.ci(j,stim,1), stats.on.ci(j,stim,2), stats.N(j,stim));
    disp(str);
end
disp(' ');

% Pentanedione OFF.
stim = 2;
neuron_i = find(stats.off.qT.off(:,stim) <= q_thresh);
neuron_off(neuron_i) = neuron_off(neuron_i) + 1;
disp(['*** 2,3-pentanedione OFF = ' num2str(length(neuron_i))]);
for i = 1:length(neuron_i)
    j = neuron_i(i);
    str = sprintf('%s (q=%.3f, ci=%d [%.2f,%.2f], N=%d)', ...
        stats.neuron{j}, stats.off.qT.off(j,stim), ...
        stats.off.ci_sign(j,stim), ...
        stats.off.ci(j,stim,1), stats.off.ci(j,stim,2), stats.N(j,stim));
    disp(str);
end
disp(' ');

% NaCl ON.
stim = 3;
neuron_i = find(stats.on.qT.on(:,stim) <= q_thresh);
neuron_on(neuron_i) = neuron_on(neuron_i) + 1;
disp(['*** NaCl ON = ' num2str(length(neuron_i))]);
for i = 1:length(neuron_i)
    j = neuron_i(i);
    str = sprintf('%s (q=%.3f, ci=%d [%.2f,%.2f], N=%d)', ...
        stats.neuron{j}, stats.on.qT.on(j,stim), ...
        stats.on.ci_sign(j,stim), ...
        stats.on.ci(j,stim,1), stats.on.ci(j,stim,2), stats.N(j,stim));
    disp(str);
end
disp(' ');

% NaCl OFF.
stim = 3;
neuron_i = find(stats.off.qT.off(:,stim) <= q_thresh);
neuron_off(neuron_i) = neuron_off(neuron_i) + 1;
disp(['*** NaCl OFF = ' num2str(length(neuron_i))]);
for i = 1:length(neuron_i)
    j = neuron_i(i);
    str = sprintf('%s (q=%.3f, ci=%d [%.2f,%.2f], N=%d)', ...
        stats.neuron{j}, stats.off.qT.off(j,stim), ...
        stats.off.ci_sign(j,stim), ...
        stats.off.ci(j,stim,1), stats.off.ci(j,stim,2), stats.N(j,stim));
    disp(str);
end
disp(' ');

% Light.
stim = 4;
neuron_i = find(stats.on.qT.on(:,stim) <= q_thresh);
neuron_on(neuron_i) = neuron_on(neuron_i) + 1;
disp(['*** light ON = ' num2str(length(neuron_i))]);
for i = 1:length(neuron_i)
    j = neuron_i(i);
    str = sprintf('%s (q=%.3f, ci= %d[%.2f,%.2f], N=%d)', ...
        stats.neuron{j}, stats.on.qT.on(j,stim), ...
        stats.on.ci_sign(j,stim), ...
        stats.on.ci(j,stim,1), stats.on.ci(j,stim,2), stats.N(j,stim));
    disp(str);
end
disp(' ');

% Show the ON neuron counts.
disp(['*** ON Neurons = ' num2str(sum(neuron_on > 0))]);
for i = 1:length(neuron_on)
    if neuron_on(i) > 0
        str = sprintf('%s = %d', stats.neuron{i}, neuron_on(i));
        disp(str);
    end
end
disp(' ');

% Show the OFF neuron counts.
disp(['*** OFF Neurons = ' num2str(sum(neuron_off > 0))]);
for i = 1:length(neuron_off)
    if neuron_off(i) > 0
        str = sprintf('%s = %d', stats.neuron{i}, neuron_off(i));
        disp(str);
    end
end
disp(' ');

% Show the ALL neuron counts.
neuron_all = neuron_on + neuron_off;
disp(['*** ALL Neurons = ' num2str(sum(neuron_all > 0))]);
for i = 1:length(neuron_all)
    if neuron_all(i) > 0
        str = sprintf('%s = %d', stats.neuron{i}, neuron_all(i));
        disp(str);
    end
end
disp(' ');
