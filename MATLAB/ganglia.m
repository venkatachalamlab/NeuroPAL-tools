% Create a list of ganglia & their neurons.
ganglion = {};

% Anterior Pharyngeal Bulb
ganglion(1).name = 'Anterior Pharyngeal Bulb';
ganglion(1).neurons = {
    'I1'
    'I2'
    'I3'
    'M3'
    'M4'
    'MC'
    'MI'
    'NSM'
    };

% Posterior Pharyngeal Bulb
ganglion(end+1).name = 'Posterior Pharyngeal Bulb';
ganglion(end).neurons = {
    'I4'
    'I5'
    'I6'
    'M1'
    'M2'
    'M5'
    };

% Anterior Ganglion
ganglion(end+1).name = 'Anterior Ganglion';
ganglion(end).neurons = {
    'BAG'
    'CEPV'
    'IL1'
    'IL1D'
    'IL1V'
    'IL2'
    'IL2D'
    'IL2V'
    'OLL'
    'OLQD'
    'OLQV'
    'RIP'
    'RME'
    'RMED'
    'RMEV'
    'URAD'
    'URAV'
    'URB'
    'URYD'
    'URYV'
    };

% Dorsal Ganglion
ganglion(end+1).name = 'Dorsal Ganglion';
ganglion(end).neurons = {
    'ALA'
    'CEPD'
    'RID'
    'URX'
    };

% Lateral Ganglion
ganglion(end+1).name = 'Lateral Ganglion';
ganglion(end).neurons = {
    'ADA'
    'ADE'
    'ADF'
    'ADL'
    'AFD'
    'AIB'
    'AIN'
    'AIZ'
    'AQR'
    'ASE'
    'ASG'
    'ASH'
    'ASI'
    'ASJ'
    'ASK'
    'AUA'
    'AVA'
    'AVB'
    'AVD'
    'AVE'
    'AVH'
    'AVJ'
    'AWA'
    'AWB'
    'AWC'
    'FLP'
    'RIA'
    'RIB'
    'RIC'
    'RIM'
    'RIV'
    'RMD'
    'RMDV'
    'RMG'
    'SAAV'
    'SIBD'
    'SMDV'
    };

% Ventral Ganglion
ganglion(end+1).name = 'Ventral Ganglion';
ganglion(end).neurons = {
    'AIA'
    'AIM'
    'AIY'
    'AVK'
    'AVL'
    'SAAD'
    'SIAD'
    'SIAV'
    'SIBV'
    'SMBD'
    'SMBV'
    'SMDD'
    'RIH'
    'RIR'
    'RIS'
    'RMDD'
    'RMF'
    'RMH'
    };

% Retrovesicular Ganglion
ganglion(end+1).name = 'Retrovesicular Ganglion';
ganglion(end).neurons = {
    'AS1'
    'AVF'
    'AVG'
    'DA1'
    'DB1'
    'DB2'
    'DD1'
    'RIF'
    'RIG'
    'SABD'
    'SABV'
    'VA1'
    'VD1'
    'VD2'
    'VB1'
    'VB2'
    };

% Midbody
ganglion(end+1).name = 'Midbody';
ganglion(end).neurons = {
    'ALM'
    'AVM'
    'BDU'
    'CAN'
    'HSN'
    'PDE'
    'PVD'
    'PVM'
    'SDQ'
    };

% Ventral Nerve Cord
ganglion(end+1).name = 'Ventral Nerve Cord';
ganglion(end).neurons = {
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
ganglion(end+1).name = 'Pre-Anal Ganglion';
ganglion(end).neurons = {
    'AS11'
    'DA8'
    'DA9'
    'DD6'
    'PDA'
    'PDB'
    'PVP'
    'PVT'
    'VA12'
    'VD12'
    'VD13'
    };

% Dorso-Rectal Ganglion
ganglion(end+1).name = 'Dorso-Rectal Ganglion';
ganglion(end).neurons = {
    'DVA'
    'DVB'
    'DVC'
    };

% Tail
ganglion(end+1).name = 'Tail';
ganglion(end).neurons = {
    'ALN'
    'LUA'
    'PHA'
    'PHB'
    'PHC'
    'PLM'
    'PLN'
    'PQR'
    'PVR'
    'PVC'
    'PVN'
    'PVQ'
    'PVW'
    };

% Save the ganglia.
save('ganglia.mat', 'ganglion');
