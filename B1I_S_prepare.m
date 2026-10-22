function result = B1I_S_prepare(cfg)
% Standalone B1I front end. No IMU values or navigation results are invented.
% No existing GPS/B1C filter or ephemeris decoder is used for B1I data.
if nargin == 0, cfg = B1I_config(); end
projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir);
if ~any(strcmp(cfg.mode, {'inspect','acquire'}))
    error('B1I:FusionUnavailable', ['This front end only supports inspect/acquire. ', ...
        'Use B1I_S_run for standalone tracking, D1/D2 decoding and positioning. ', ...
        'IMU deep fusion is not implemented.']);
end
validateattributes(cfg.sampleFreq, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(cfg.startSeconds, {'numeric'}, {'scalar','real','finite','nonnegative'});
assert(strcmp(cfg.iqFormat, 'fc32'), 'B1I:Format', 'This reader expects fc32 IQ.');
iqInfo = dir(cfg.iqFile);
assert(numel(iqInfo)==1 && ~iqInfo.isdir, 'B1I:MissingIQ', 'IQ file not found.');
result = struct('stage','data inspection', 'fusionPerformed',false, 'config',cfg);
result.iq.bytes = iqInfo.bytes;
result.iq.completeComplexSamples = floor(iqInfo.bytes/8);
result.iq.trailingBytes = mod(iqInfo.bytes,8);
result.iq.seconds = result.iq.completeComplexSamples/cfg.sampleFreq;
fprintf('B1I IQ duration: %.6f s at %.0f samples/s (fc32).\n', ...
    result.iq.seconds, cfg.sampleFreq);
if result.iq.trailingBytes ~= 0
    warning('B1I:TruncatedIQ', ['IQ ends with %d incomplete-sample bytes. ', ...
        'They are ignored for acquisition; check the original recording/copy.'], ...
        result.iq.trailingBytes);
end
result.imu = BDS.B1I.inspect_imu(cfg.imuFile);
fprintf('IMU inspection: %s\n', result.imu.message);
if ~cfg.sync.confirmed
    fprintf('IMU/IQ synchronization is unconfirmed. Fusion is disabled.\n');
end
if strcmp(cfg.mode, 'acquire')
    result.stage = 'B1I acquisition only';
    result.acquisition = BDS.B1I.acquisition(cfg);
    disp(result.acquisition);
    fprintf('%d candidate(s); acquisition is NOT tracking or a position solution.\n', ...
        nnz(result.acquisition.Candidate));
    if cfg.showPlot
        figure('Name','B1I acquisition candidates - no fusion');
        bar(result.acquisition.PRN, result.acquisition.PeakPowerRatio);
        hold on;
        plot([0,64], cfg.peakRatioThreshold*[1,1], 'r--');
        xlabel('BDS PRN'); ylabel('Peak / second peak (power)'); grid on;
        title('B1I acquisition only; candidate threshold is heuristic');
    end
end
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
resultFile = [tempname(cfg.outputDir), '.mat'];
save(resultFile, 'result', 'cfg');
fprintf('Saved inspection/acquisition result: %s\n', resultFile);
end
