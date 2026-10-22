function B1I_selftest()
% Run in MATLAB to check the MATLAB implementation. No hardware data needed.
% This checks acquisition, not tracking, navigation decoding or fusion.
cfg = B1I_config();
cfg.prns = [19,2];
cfg.noncoherentPeriods = 8;
cfg.dopplerMaxHz = 2500;
cfg.dopplerStepHz = 250;
cfg.ifHz = 0;
cfg.iqOrder = 'IQ';
cfg.byteOrder = 'ieee-le';
cfg.startSeconds = 0;
cfg.showPlot = false;
cfg.iqFile = [tempname, '.dat'];
cleanup = onCleanup(@() remove_fixture(cfg.iqFile)); %#ok<NASGU>
fs = cfg.sampleFreq;
n = fs/1000;
samples = 0:(2*n*cfg.noncoherentPeriods-1);
delay = 1234;
frequency = 1750;
primary = BDS.B1I.codeGene(19);
chips = floor((samples-delay)*2.046e6/fs);
period = floor((samples-delay)/n);
% Alternate every primary-code period to exercise boundary sign changes.
signal = primary(mod(chips,2046)+1) .* (-1).^period;
signal = signal .* exp(2i*pi*frequency*samples/fs);
fid = fopen(cfg.iqFile,'wb','ieee-le');
assert(fid>=0);
closeFixture = onCleanup(@() fclose(fid));
fwrite(fid,[real(signal);imag(signal)],'single');
clear closeFixture;
readback = BDS.B1I.read_iq(cfg.iqFile,0,numel(signal),'IQ','ieee-le');
assert(max(abs(signal-readback))<1e-6, 'fc32 readback failed.');
shifted = BDS.B1I.read_iq(cfg.iqFile,17,10,'IQ','ieee-le');
assert(max(abs(shifted-signal(18:27)))<1e-6, 'Byte offset failed.');
swapped = BDS.B1I.read_iq(cfg.iqFile,0,10,'QI','ieee-le');
assert(max(abs(swapped-(imag(signal(1:10))+1i*real(signal(1:10)))))<1e-6);
assert_error(@() BDS.B1I.read_iq(cfg.iqFile,0,numel(signal)+1,'IQ','ieee-le'), ...
    'B1I:ShortIQ');
result = BDS.B1I.acquisition(cfg);
row = result(result.PRN==19,:);
assert(abs(row.CodeDelaySamples-delay)<=1, 'Incorrect code delay.');
assert(row.DopplerHz==frequency, 'Incorrect Doppler.');
assert(row.Candidate, 'Known signal not detected.');
for prn = 1:63
    primary = BDS.B1I.codeGene(prn);
    assert(numel(primary)==2046 && all(abs(primary)==1));
end
assert_error(@() BDS.B1I.codeGene(0), '');
assert_error(@() BDS.B1I.codeGene(64), '');
cfg.mode = 'fusion';
assert_error(@() B1I_S_prepare(cfg), 'B1I:FusionUnavailable');
fprintf('B1I MATLAB self-test passed (front end only).\n');
end

function assert_error(action, expected)
try
    action();
catch err
    if ~isempty(expected), assert(strcmp(err.identifier,expected)); end
    return;
end
error('B1I:TestExpectedError','Expected an error, but the operation succeeded.');
end

function remove_fixture(filename)
if isfile(filename), delete(filename); end
end
