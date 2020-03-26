% Create a list of ganglia & their neurons.
ganglion_LR = {};

% Anterior Pharyngeal Bulb
ganglion_LR(1).name = 'Anterior Pharyngeal Bulb';
ganglion_LR(1).neurons = {
    'I1L'
    'I1R'
    'I2L'
    'I2R'
    'I3'
    'M3L'
    'M3R'
    'M4'
    'MCL'
    'MCR'
    'MI'
    'NSML'
    'NSMR'
    };

% Posterior Pharyngeal Bulb
ganglion_LR(end+1).name = 'Posterior Pharyngeal Bulb';
ganglion_LR(end).neurons = {
    'I4'
    'I5'
    'I6'
    'M1'
    'M2L'
    'M2R'
    'M5'
    };

% Anterior Ganglion (Left)
ganglion_LR(end+1).name = 'Anterior Ganglion (Left)';
ganglion_LR(end).neurons = {
    'BAGL'
    'CEMVL'
    'CEPVL'
    'IL1L'
    'IL1DL'
    'IL1VL'
    'IL2L'
    'IL2DL'
    'IL2VL'
    'MCML'
    'OLLL'
    'OLQDL'
    'OLQVL'
    'RIPL'
    'RMEL'
    'RMED'
    'RMEV'
    'URADL'
    'URAVL'
    'URBL'
    'URYDL'
    'URYVL'
    };

% Anterior Ganglion (Right)
ganglion_LR(end+1).name = 'Anterior Ganglion (Right)';
ganglion_LR(end).neurons = {
    'BAGR'
    'CEMVR'
    'CEPVR'
    'IL1R'
    'IL1DR'
    'IL1VR'
    'IL2R'
    'IL2DR'
    'IL2VR'
    'MCMR'
    'OLLR'
    'OLQDR'
    'OLQVR'
    'RIPR'
    'RMER'
    'RMED'
    'RMEV'
    'URADR'
    'URAVR'
    'URBR'
    'URYDR'
    'URYVR'
    };

% Dorsal Ganglion
ganglion_LR(end+1).name = 'Dorsal Ganglion';
ganglion_LR(end).neurons = {
    'ALA'
    'CEMDL'
    'CEMDR'
    'CEPDL'
    'CEPDR'
    'RID'
    'URXL'
    'URXR'
    };

% Lateral Ganglion (Left)
ganglion_LR(end+1).name = 'Lateral Ganglion (Left)';
ganglion_LR(end).neurons = {
    'ADAL'
    'ADEL'
    'ADFL'
    'ADLL'
    'AFDL'
    'AIBL'
    'AINL'
    'AIZL'
    'ASEL'
    'ASGL'
    'ASHL'
    'ASIL'
    'ASJL'
    'ASKL'
    'AUAL'
    'AVAL'
    'AVBL'
    'AVDL'
    'AVEL'
    'AVHL'
    'AVJL'
    'AWAL'
    'AWBL'
    'AWCL'
    'AWCL(OFF)'
    'AWCL(ON)'
    'FLPL'
    'RIAL'
    'RIBL'
    'RICL'
    'RIML'
    'RIVL'
    'RMDL'
    'RMDVL'
    'RMGL'
    'SAAVL'
    'SIBDL'
    'SMDVL'
    };

% Lateral Ganglion (Right)
ganglion_LR(end+1).name = 'Lateral Ganglion (Right)';
ganglion_LR(end).neurons = {
    'ADAR'
    'ADER'
    'ADFR'
    'ADLR'
    'AFDR'
    'AIBR'
    'AINR'
    'AIZR'
    'AQR'
    'ASER'
    'ASGR'
    'ASHR'
    'ASIR'
    'ASJR'
    'ASKR'
    'AUAR'
    'AVAR'
    'AVBR'
    'AVDR'
    'AVER'
    'AVHR'
    'AVJR'
    'AWAR'
    'AWBR'
    'AWCR'
    'AWCR(OFF)'
    'AWCR(ON)'
    'FLPR'
    'RIAR'
    'RIBR'
    'RICR'
    'RIMR'
    'RIVR'
    'RMDR'
    'RMDVR'
    'RMGR'
    'SAAVR'
    'SIBDR'
    'SMDVR'
    };

% Ventral Ganglion
ganglion_LR(end+1).name = 'Ventral Ganglion';
ganglion_LR(end).neurons = {
    'AIAL'
    'AIAR'
    'AIML'
    'AIMR'
    'AIYL'
    'AIYR'
    'AVKL'
    'AVKR'
    'AVL'
    'SAADL'
    'SAADR'
    'SIADL'
    'SIADR'
    'SIAVL'
    'SIAVR'
    'SIBVL'
    'SIBVR'
    'SMBDL'
    'SMBDR'
    'SMBVL'
    'SMBVR'
    'SMDDL'
    'SMDDR'
    'RIH'
    'RIR'
    'RIS'
    'RMDDL'
    'RMDDR'
    'RMFL'
    'RMFR'
    'RMHL'
    'RMHR'
    };

% Retrovesicular Ganglion
ganglion_LR(end+1).name = 'Retro-Vesicular Ganglion';
ganglion_LR(end).neurons = {
    'AS1'
    'AVFL'
    'AVFR'
    'AVG'
    'DA1'
    'DB1'
    'DB2'
    'DD1'
    'RIFL'
    'RIFR'
    'RIGL'
    'RIGR'
    'SABD'
    'SABVL'
    'SABVR'
    'VA1'
    'VD1'
    'VD2'
    'VB1'
    'VB2'
    };

% Anterior Midbody
ganglion_LR(end+1).name = 'Anterior Midbody';
ganglion_LR(end).neurons = {
    'ALML'
    'ALMR'
    'AVM'
    'BDUL'
    'BDUR'
    'SDQR'
    };

% Central Midbody
ganglion_LR(end+1).name = 'Central Midbody';
ganglion_LR(end).neurons = {
    'CANL'
    'CANR'
    };

% Posterior Midbody
ganglion_LR(end+1).name = 'Posterior Midbody';
ganglion_LR(end).neurons = {
    'PDEL'
    'PDER'
    'PVDL'
    'PVDR'
    'PVM'
    'SDQL'
    };

% Ventral Nerve Cord
ganglion_LR(end+1).name = 'Ventral Nerve Cord';
ganglion_LR(end).neurons = {
    'AS2'
    'AS3'
    'AS4'
    'AS5'
    'AS6'
    'AS7'
    'AS8'
    'AS9'
    'AS10'
    'CA1'
    'CA2'
    'CA3'
    'CA4'
    'CA5'
    'CA6'
    'CA7'
    'CA8'
    'CA9'
    'CP1'
    'CP2'
    'CP3'
    'CP4'
    'CP5'
    'CP6'
    'CP7'
    'CP8'
    'DA2'
    'DA3'
    'DA4'
    'DA5'
    'DA6'
    'DA7'
    'DB3'
    'DB4'
    'DB5'
    'DB6'
    'DB7'
    'DD2'
    'DD3'
    'DD4'
    'DD5'
    'VA2'
    'VA3'
    'VA4'
    'VA5'
    'VA6'
    'VA7'
    'VA8'
    'VA9'
    'VA10'
    'VA11'
    'VB3'
    'VB4'
    'VB5'
    'VB6'
    'VB7'
    'VB8'
    'VB9'
    'VB10'
    'VB11'
    'VD3'
    'VD4'
    'VD5'
    'VD6'
    'VD7'
    'VD8'
    'VD9'
    'VD10'
    'VD11'
    };

% Pre-Anal Ganglion
ganglion_LR(end+1).name = 'Pre-Anal Ganglion';
ganglion_LR(end).neurons = {
    'AS11'
    'CP9'
    'DA8'
    'DA9'
    'DD6'
    'HOA'
    'HOB'
    'PDA'
    'PDB'
    'PDC'
    'PGA'
    'PVP'
    'PVS'
    'PVT'
    'PVU'
    'PVV'
    'PVX'
    'PVY'
    'PVZ'
    'VA12'
    'VD12'
    'VD13'
    };

% Dorso-Rectal Ganglion
ganglion_LR(end+1).name = 'Dorso-Rectal Ganglion';
ganglion_LR(end).neurons = {
    'DVA'
    'DVB'
    'DVC'
    'DVE'
    'DVF'
    'DX1'
    'DX2'
    'DX3'
    'DX4'
    'EF1'
    'EF2'
    'EF3'
    'EF4'
    };

% Rays (Left)
ganglion_LR(end+1).name = 'Rays (Left)';
ganglion_LR(end).neurons = {
    'R1AL'
    'R1BL'
    'R2AL'
    'R2BL'
    'R3AL'
    'R3BL'
    'R4AL'
    'R4BL'
    'R5AL'
    'R5BL'
    'R6AL'
    'R6BL'
    'R7AL'
    'R7BL'
    'R8AL'
    'R8BL'
    'R9AL'
    'R9BL'
    };

% Rays (Right)
ganglion_LR(end+1).name = 'Rays (Right)';
ganglion_LR(end).neurons = {
    'R1AR'
    'R1BR'
    'R2AR'
    'R2BR'
    'R3AR'
    'R3BR'
    'R4AR'
    'R4BR'
    'R5AR'
    'R5BR'
    'R6AR'
    'R6BR'
    'R7AR'
    'R7BR'
    'R8AR'
    'R8BR'
    'R9AR'
    'R9BR'
    };

% Cloacal Ganglion
ganglion_LR(end+1).name = 'Cloacal Ganglion (Left)';
ganglion_LR(end).neurons = {
    'PCAL'
    'PCBL'
    'PCCL'
    'SPCL'
    'SPDL'
    'SPVL'
    };

% Cloacal Ganglion
ganglion_LR(end+1).name = 'Cloacal Ganglion (Right)';
ganglion_LR(end).neurons = {
    'PCAR'
    'PCBR'
    'PCCR'
    'SPCR'
    'SPDR'
    'SPVR'
    };

% Tail (Left)
ganglion_LR(end+1).name = 'Lumbar Ganglion (Left)';
ganglion_LR(end).neurons = {
    'ALNL'
    'LUAL'
    'PHAL'
    'PHBL'
    'PHCL'
    'PHDL'
    'PLML'
    'PLNL'
    'PQR'
    'PVCL'
    'PVNL'
    'PVQL'
    'PVWL'
    };

% Tail (Right)
ganglion_LR(end+1).name = 'Lumbar Ganglion (Right)';
ganglion_LR(end).neurons = {
    'ALNR'
    'LUAR'
    'PHAR'
    'PHBR'
    'PHCR'
    'PHDR'
    'PLMR'
    'PLNR'
    'PVCR'
    'PVNR'
    'PVQR'
    'PVR'
    'PVWR'
    };

% Save the ganglia.
save('ganglia_male_LR.mat', 'ganglion_LR');
