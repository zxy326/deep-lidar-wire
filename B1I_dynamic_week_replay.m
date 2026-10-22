function report=B1I_dynamic_week_replay(dataDir)
% End-to-end delayed sixth satellite and week rollover using generated IQ.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
m=jsondecode(fileread(fullfile(dataDir,'manifest.json')));
assert(m.latePrnSeconds>24 && m.receiverBdtStart+m.config.processSeconds>604800);
cfg=merge(B1I_deep_config(),m.config); cfg.showPlot=false;
cfg.prns=cfg.prns(:)'; cfg.deep.initialAttitudeDeg=cfg.deep.initialAttitudeDeg(:)';
cfg.deep.armBodyMeters=cfg.deep.armBodyMeters(:)';
r=BDS.B1I.deep_receive(cfg);
assert(r.timeSeconds(end)>=cfg.startSeconds+cfg.processSeconds-cfg.deep.dt/2,'Replay stopped early.');
j=find(r.prns==m.latePrn); assert(isscalar(j),'New satellite was not discovered.');
assert(r.channelAddedTimeSeconds(j)>=m.latePrnSeconds,'Absent satellite admitted early.');
used=r.diagnostics.used(:,j);
assert(any(used),'New satellite never contributed a validated measurement.');
first=r.timeSeconds(find(used,1));
assert(first-r.channelAddedTimeSeconds(j)>12,'New satellite skipped fresh-navigation quarantine.');
assert(any(r.feedbackApplied(:,j)),'New satellite never received closed-loop feedback.');
assert(any(r.bdtWeek==m.bdtWeek) && any(r.bdtWeek==m.bdtWeek+1));
continuous=(r.bdtWeek-m.bdtWeek)*604800+r.receiverBdtSow;
assert(max(abs(diff(continuous)-diff(r.timeSeconds)))<1e-8,'Receiver time jumped across week.');
assert(numel(unique(r.prns))==numel(r.prns));
truth=readmatrix(fullfile(dataDir,'truth.csv'),'NumHeaderLines',1);
ref=interp1(truth(:,1),truth(:,8:13),r.timeSeconds,'linear');
report=struct('passed',true,'runDirectory',r.runDir,'newPrn',m.latePrn, ...
    'discoveredSeconds',r.channelAddedTimeSeconds(j),'firstUsedSeconds',first, ...
    'newPrnUpdateEpochs',nnz(used),'epochs',numel(r.timeSeconds), ...
    'heightRmseMeters',sqrt(mean((r.positionLla(:,3)-ref(:,3)).^2)), ...
    'velocityRmseMps',sqrt(mean(sum((r.velocityNed-ref(:,4:6)).^2,2))), ...
    'note','Functional synthetic D1 acceptance; not real data or D2 rollover acceptance.');
save(fullfile(r.runDir,'dynamic_week_report.mat'),'report'); disp(report);
end
function cfg=merge(cfg,x)
for name=fieldnames(x)'
    n=name{1}; if isstruct(x.(n))&&isfield(cfg,n), cfg.(n)=merge(cfg.(n),x.(n)); else, cfg.(n)=x.(n); end
end
end
