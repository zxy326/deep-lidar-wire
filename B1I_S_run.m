%% Independent BeiDou B1I receiver: acquisition -> tracking -> navigation.
% Original GL1CA_S_deep.m and all GPS receiver classes are unchanged.
projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir);
b1iResult = BDS.B1I.receive(B1I_config());
