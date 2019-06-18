% Create a list of stimuli and assign them colors.
i = 1;
stimulus = [];

% 2-Butanone.
stimulus(i).name = '10^{-4} 2-butanone';
stimulus(i).keyword = 'butanone';
stimulus(i).RGB = [0,1,1];
i = i + 1;

% 2,3-Pentanedione.
stimulus(i).name = '10^{-4} 2,3-pentanedione';
stimulus(i).keyword = 'pentanedione';
stimulus(i).RGB = [1,0,1];
i = i + 1;

% NaCl.
stimulus(i).name = '200uM NaCl';
stimulus(i).keyword = 'nacl';
stimulus(i).RGB = [1,1,0];
i = i + 1;

% Save the stimuli.
save('stimuli.mat', 'stimulus');
