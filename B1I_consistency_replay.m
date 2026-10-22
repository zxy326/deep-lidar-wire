function report=B1I_consistency_replay(dataDir,options)
% One retained fixture: full default replay, short feedback-disabled smoke test.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
manifest=jsondecode(fileread(fullfile(dataDir,'manifest.json')));
cfg=merge(B1I_deep_config(),manifest.config); cfg.showPlot=false;
cfg.prns=cfg.prns(:)'; cfg.deep.initialAttitudeDeg=cfg.deep.initialAttitudeDeg(:)';
cfg.deep.armBodyMeters=cfg.deep.armBodyMeters(:)';
cfg.deep.physicalRangeRate=true; cfg.deep.sourceCodeVariance=false;
cfg.deep.feedbackGain=1;
if nargin>1, cfg.deep=merge(cfg.deep,options); end
current=BDS.B1I.deep_receive(cfg);
assert(current.closedLoop && current.filterPhysicalRangeRate);
assert(all(current.feedbackApplied(:)==current.feedbackEligible(:)));
assert(~any(current.feedbackNonzero(:)&~current.feedbackApplied(:)));
assert(current.timeSeconds(end)>=cfg.startSeconds+cfg.processSeconds-cfg.deep.dt/2);
lag=current.diagnostics.channelBoundaryLagSeconds;
assert(all(lag(:)>=-1e-9 & lag(:)<.002));
summary=jsondecode(fileread(fullfile(dataDir,'physical_rate_report.json')));
old=BDS.B1I.load_result(fullfile(summary.C.runDirectory,'deep_result.mat'));
truth=readmatrix(fullfile(dataDir,'truth.csv'),'NumHeaderLines',1);
report=struct('previous',metrics(old,truth),'current',metrics(current,truth), ...
    'runDirectory',current.runDir,'maxChannelLagSeconds',max(lag(:)), ...
    'note','Single-fixture regression; no parameter tuning; not a multi-scenario accuracy claim.');
report.initializationInfo=current.initializationInfo;
report.navigationRefreshCount=current.navigationRefreshCount;
off=cfg; off.deep.feedbackGain=0;
off.processSeconds=current.initialTimeSeconds-off.startSeconds+1.25;
control=BDS.B1I.deep_receive(off);
assert(~control.closedLoop && ~any(control.feedbackApplied(:)) && ...
    ~any(control.feedbackNonzero(:)) && any(control.feedbackEligible(:)));
assert(~any(control.codeCorrectionChips(:)) && ~any(control.carrierCorrectionHz(:)));
assert(control.timeSeconds(end)>=off.startSeconds+off.processSeconds-off.deep.dt/2);
report.disabledFeedbackRunDirectory=control.runDir;
report.disabledFeedbackEpochs=numel(control.timeSeconds);
file=fullfile(current.runDir,'consistency_replay.json');
fid=fopen(file,'w'); assert(fid>=0); clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
disp(report); fprintf('Consistency replay report: %s\n',file);
end
function m=metrics(r,truth)
ref=interp1(truth(:,1),truth(:,8:13),r.timeSeconds,'linear');
dh=r.positionLla(:,3)-ref(:,3); dv=r.velocityNed(:,3)-ref(:,6);
assert(all(isfinite([dh;dv])));
dp=zeros(numel(dh),3);
for k=1:numel(dh), dp(k,:)=lla2ecef(r.positionLla(k,:))-lla2ecef(ref(k,1:3)); end
velocityError=r.velocityNed-ref(:,4:6);
m=struct('heightRmseMeters',sqrt(mean(dh.^2)), ...
 'maxHeightMeters',max(abs(dh)),'downVelocityRmseMps',sqrt(mean(dv.^2)), ...
 'positionRmseMeters',sqrt(mean(sum(dp.^2,2))), ...
 'velocityRmseMps',sqrt(mean(sum(velocityError.^2,2))), ...
 'epochs',numel(dh),'updates',nnz(r.measurementUsed),'stop',r.stopReason);
if isfield(r.diagnostics,'beforeCorrection')
 correction=r.diagnostics.afterCorrection(:,3)-r.diagnostics.beforeCorrection(:,3);
 m.firstHeightCorrectionMeters=correction(1);
 m.firstSecondMaxHeightCorrectionMeters=max(abs(correction(r.timeSeconds<=r.initialTimeSeconds+1)));
end
end
function cfg=merge(cfg,x)
for name=fieldnames(x)'
 n=name{1}; if isstruct(x.(n))&&isfield(cfg,n), cfg.(n)=merge(cfg.(n),x.(n)); else, cfg.(n)=x.(n); end
end
end
