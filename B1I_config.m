function cfg = B1I_config()
% Settings for the user's B1I recording. ASCII source avoids legacy encoding
% issues. The directory below is D:\deep followed by the Chinese word zuhe.
cfg.dataDir = fullfile('D:\', ['deep', char([32452,21512])]);
cfg.iqFile = fullfile(cfg.dataDir, 'bds_b1i_20260921_5Msps_fc32.dat');
cfg.imuFile = fullfile(cfg.dataDir, 'deta40_raw.bag');
cfg.sampleFreq = 5e6;             % Hz; inferred from filename, VERIFY.
cfg.iqFormat = 'fc32';            % I,Q interleaved IEEE single precision.
cfg.byteOrder = 'ieee-le';        % Little endian; VERIFY with recorder.
cfg.iqOrder = 'IQ';              % Change to QI if the recorder writes Q,I.
cfg.ifHz = 0;                    % Signal IF in complex IQ; VERIFY, not RF Hz.
cfg.startSeconds = 0;
cfg.mode = 'acquire';            % 'inspect' or 'acquire'; fusion unavailable.
cfg.prns = 1:63;
cfg.dopplerMaxHz = 10000;
cfg.dopplerStepHz = 250;
cfg.noncoherentPeriods = 20;     % 20 separate 1-ms integrals in 40-ms input.
cfg.peakRatioThreshold = 2.5;    % Heuristic POWER ratio; not a lock guarantee.
cfg.showPlot = true;
cfg.outputDir = fullfile(fileparts(mfilename('fullpath')), '~temp', 'b1i');

% Future fusion prerequisites. No automatic alignment of the two first rows.
cfg.sync.confirmed = false;
cfg.sync.iqStartUtc = [];        % Exact UTC [Y M D h m s] for first IQ sample.
cfg.sync.imuToIqOffsetSeconds = NaN; % Must be measured, not guessed.
cfg.sync.imuAxesToBody = [];     % 3x3 rotation, established from mounting.
% Standalone receiver settings (B1I_S_run); IMU deep integration is separate.
cfg.processSeconds = 100;
cfg.maxChannels = 12;
cfg.pllBandwidthHz = 18;
cfg.dllBandwidthHz = 1;
cfg.codeReplica='rectangular'; % Real-data front-end response must be established separately.
cfg.trackCacheSeconds = 0.25;
cfg.maxBchCorrections = 0;       % Strict parity by default (no silent repairs).
cfg.navStepSeconds = 1;
cfg.maxEphemerisAgeSeconds = 7200;
cfg.maxResidualMeters = 150;
cfg.maxGdop = 30;
cfg.minElevationDeg = 5;
cfg.initialEcef = [0,0,0];       % Earth-center startup; no guessed location.
cfg.positionAtmosphere = true;  % Dry troposphere only; no ionosphere correction.
cfg.acquisitionScanOffsets = [0,5,10,20,30]; % Relative to startSeconds.
cfg.probeSeconds = 3;           % Require sustained tracking before long run.
end
