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
    'CEPVL'
    'IL1L'
    'IL1DL'
    'IL1VL'
    'IL2L'
    'IL2DL'
    'IL2VL'
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
    'CEPVR'
    'IL1R'
    'IL1DR'
    'IL1VR'
    'IL2R'
    'IL2DR'
    'IL2VR'
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
ganglion_LR(end+1).name = 'Retrovesicular Ganglion';
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

% Midbody (Left)
ganglion_LR(end+1).name = 'Midbody (Left)';
ganglion_LR(end).neurons = {
    'ALML'
    'BDUL'
    'CANL'
    'HSNL'
    'PDEL'
    'PVDL'
    'PVM'
    'SDQL'
    };

% Midbody (Right)
ganglion_LR(end+1).name = 'Midbody (Right)';
ganglion_LR(end).neurons = {
    'ALMR'
    'AVM'
    'BDUR'
    'CANR'
    'HSNR'
    'PDER'
    'PVDR'
    'SDQR'
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
    'VC1'
    'VC2'
    'VC3'
    'VC4'
    'VC5'
    'VC6'
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
    'DA8'
    'DA9'
    'DD6'
    'PDA'
    'PDB'
    'PVPL'
    'PVPR'
    'PVT'
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
    };

% Tail (Left)
ganglion_LR(end+1).name = 'Tail (Left)';
ganglion_LR(end).neurons = {
    'ALNL'
    'LUAL'
    'PHAL'
    'PHBL'
    'PHCL'
    'PLML'
    'PLNL'
    'PQR'
    'PVCL'
    'PVNL'
    'PVQL'
    'PVWL'
    };

% Tail (Right)
ganglion_LR(end+1).name = 'Tail (Right)';
ganglion_LR(end).neurons = {
    'ALNR'
    'LUAR'
    'PHAR'
    'PHBR'
    'PHCR'
    'PLMR'
    'PLNR'
    'PVCR'
    'PVNR'
    'PVQR'
    'PVR'
    'PVWR'
    };

% Save the ganglia.
save('ganglia_LR.mat', 'ganglion_LR');
