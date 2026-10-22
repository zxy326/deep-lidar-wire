function report=B1I_frozen_ablation(dataDir,baselineFile,mode)
% Freeze saved baseline, then change exactly one named factor on the same IQ/IMU.
% mode: baseline, range, rate, rate_late (fixed 1-s startup intervention).
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
assert(any(strcmp(mode,{'baseline','range','rate','rate_late'})));
B1I_source_quality_selftest();
B=BDS.B1I.load_result(baselineFile);
[a,b,c]=BDS.B1I.quality_options(B.config.deep); assert(~a&&~b&&~c);
cfg=B.config; cfg.showPlot=false;
cfg.deep.sourceQualityModel=false; cfg.deep.sourceRangeWeight=strcmp(mode,'range');
cfg.deep.sourceRateWeight=any(strcmp(mode,{'rate','rate_late'})); cfg.deep.sourceCn0Selection=false;
if isfield(cfg.deep,'sourceWeightStartIqSeconds'), cfg.deep=rmfield(cfg.deep,'sourceWeightStartIqSeconds'); end
if strcmp(mode,'rate_late'), cfg.deep.sourceWeightStartIqSeconds=B.initialTimeSeconds+1; end
assert(~cfg.deep.sourceCodeVariance,'Ablation baseline must use fixed code variance.');
parent=fullfile(dataDir,'frozen_ablation'); if ~isfolder(parent), mkdir(parent); end
out=tempname(parent); mkdir(out);
copyfile(baselineFile,fullfile(out,'frozen_baseline.mat'));
copyfile(fullfile(dataDir,'manifest.json'),fullfile(out,'frozen_manifest.json'));
save(fullfile(out,'trial_config.mat'),'cfg','mode');
diary(fullfile(out,'trial.log')); cleanup=onCleanup(@()diary('off'));
C=BDS.B1I.deep_receive(cfg);
assert(abs(B.initialTimeSeconds-C.initialTimeSeconds)<1e-9);
assert(max(abs(B.diagnostics.epochStartState(1,:)-C.diagnostics.epochStartState(1,:)))<1e-9);
truth=readmatrix(fullfile(dataDir,'truth.csv'),'NumHeaderLines',1);
report=struct('mode',mode,'baselineFile',baselineFile,'outputDirectory',out, ...
    'baseline',metrics(B,truth),'candidate',metrics(C,truth), ...
    'sameInitialState',true,'startupDefinition','First 1 second after fusion start', ...
    'stableDefinition','At least 2 seconds after fusion start');
report.nonRegression=report.candidate.completed && ...
    report.candidate.positionRmse<=report.baseline.positionRmse*1.05 && ...
    report.candidate.velocityRmse<=report.baseline.velocityRmse*1.05 && ...
    report.candidate.heightPeak<=report.baseline.heightPeak*1.05 && ...
    report.candidate.positionPeak<=report.baseline.positionPeak*1.05 && ...
    report.candidate.unlockedFraction<=report.baseline.unlockedFraction+.001 && ...
    report.candidate.updates>=report.baseline.updates;
report.criterion='Predeclared engineering screen: complete, updates not fewer, pos/vel RMSE and height/pos peak <=105% baseline; unlock fraction <=baseline+0.001. Not universal acceptance.';
fid=fopen(fullfile(out,'report.json'),'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fprintf('%s\n',jsonencode(report));
fprintf('ABLATION REPORT: %s\n',out);
end
function m=metrics(r,truth)
t=r.timeSeconds; ref=interp1(truth(:,1),truth(:,8:13),t);
h=r.positionLla(:,3)-ref(:,3); d=r.velocityNed(:,3)-ref(:,6);
dp=zeros(numel(t),3);
for k=1:numel(t), dp(k,:)=lla2ecef(r.positionLla(k,:))-lla2ecef(ref(k,1:3)); end
dv=r.velocityNed-ref(:,4:6); pe=vecnorm(dp,2,2); ve=vecnorm(dv,2,2);
m=struct('runDirectory',r.runDir,'epochs',numel(t),'updates',nnz(r.measurementUsed), ...
    'completed',t(end)>=r.config.startSeconds+r.config.processSeconds-.005, ...
    'stopReason',r.stopReason,'heightRmse',sqrt(mean(h.^2)),'heightPeak',max(abs(h)), ...
    'downRmse',sqrt(mean(d.^2)),'positionRmse',sqrt(mean(pe.^2)), ...
    'positionPeak',max(pe),'velocityRmse',sqrt(mean(ve.^2)), ...
    'unlockedFraction',nnz(r.diagnostics.reason==4)/numel(r.diagnostics.reason));
if isfield(r,'cn0DbHz'), m.minimumCn0=min(r.cn0DbHz(:)); end
for item={'startup','stable'}
 name=item{1}; ix=t<=r.initialTimeSeconds+1;
 if strcmp(name,'stable'), ix=t>=r.initialTimeSeconds+2; end
 m.(name)=struct('epochs',nnz(ix),'heightRmse',sqrt(mean(h(ix).^2)), ...
    'heightMean',mean(h(ix)),'heightPeak',max(abs(h(ix))), ...
    'downRmse',sqrt(mean(d(ix).^2)),'downMean',mean(d(ix)), ...
    'positionRmse',sqrt(mean(pe(ix).^2)),'positionPeak',max(pe(ix)), ...
    'velocityRmse',sqrt(mean(ve(ix).^2)));
end
% Startup before/after filter corrections: diagnose jump without changing initialization.
z=r.diagnostics; ix=t<=r.initialTimeSeconds+1;
if isfield(z,'beforeCorrection')
 m.startup.maxHeightCorrection=max(abs(z.afterCorrection(ix,3)-z.beforeCorrection(ix,3)));
 m.startup.maxDownCorrection=max(abs(z.afterCorrection(ix,6)-z.beforeCorrection(ix,6)));
end
end
