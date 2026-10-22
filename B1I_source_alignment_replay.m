function report=B1I_source_alignment_replay(dataDir,baselineFile)
% Full IQ replay of source-interface alignment, without tuning parameters.
% Optional baselineFile is a saved pre-refactor deep_result.mat.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
m=jsondecode(fileread(fullfile(dataDir,'manifest.json')));
cfg=merge(B1I_deep_config(),m.config); cfg.showPlot=false;
cfg.prns=cfg.prns(:)'; cfg.deep.initialAttitudeDeg=cfg.deep.initialAttitudeDeg(:)';
cfg.deep.armBodyMeters=cfg.deep.armBodyMeters(:)';
r=BDS.B1I.deep_receive(cfg); N=numel(r.timeSeconds);
assert(N>0 && r.timeSeconds(end)>=cfg.startSeconds+cfg.processSeconds-cfg.deep.dt/2);
assert(r.resultSchemaVersion==3 && all(size(r.satnav)==[N,8]));
assert(all(size(r.attitudeSigmaDeg)==[N,3]));
assert(all(isfinite(r.antennaPositionLla(:))) && all(isfinite(r.attitudeSigmaDeg(:))));
assert(any(all(isfinite(r.satnav),2)),'Same-loop GNSS solution was never available.');
assert(all(r.feedbackApplied(:)==r.feedbackEligible(:)));
if norm(cfg.deep.armBodyMeters)==0
    assert(isequal(r.antennaPositionLla,r.positionLla)&&isequal(r.antennaVelocityNed,r.velocityNed));
end
truth=readmatrix(fullfile(dataDir,'truth.csv'),'NumHeaderLines',1);
ref=interp1(truth(:,1),truth(:,8:13),r.timeSeconds,'linear');
dh=r.positionLla(:,3)-ref(:,3); dv=r.velocityNed-ref(:,4:6);
report=struct('runDirectory',r.runDir,'epochs',N,'updates',nnz(r.measurementUsed), ...
 'sameLoopFixes',nnz(all(isfinite(r.satnav),2)), ...
 'heightRmseMeters',sqrt(mean(dh.^2)),'velocityRmseMps',sqrt(mean(sum(dv.^2,2))), ...
 'note','Source-interface refactor; same-loop GNSS diagnostics are NOT an independent control.');
if nargin>1
    old=BDS.B1I.load_result(baselineFile);
    assert(numel(old.timeSeconds)==N && max(abs(old.timeSeconds-r.timeSeconds))<1e-8);
    dp=zeros(N,3);
    for k=1:N, dp(k,:)=lla2ecef(r.positionLla(k,:))-lla2ecef(old.positionLla(k,:)); end
    report.baselineFile=baselineFile;
    report.maxPositionChangeMeters=max(vecnorm(dp,2,2));
    report.maxVelocityChangeMps=max(vecnorm(r.velocityNed-old.velocityNed,2,2));
    assert(report.maxPositionChangeMeters<.001 && report.maxVelocityChangeMps<1e-5, ...
        'B1I:AlignmentRegression','Interface refactor changed navigation beyond numerical tolerance.');
end
file=fullfile(r.runDir,'source_alignment_report.json');
fid=fopen(file,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true)); disp(report); fprintf('Report: %s\n',file);
end
function cfg=merge(cfg,x)
for name=fieldnames(x)'
 n=name{1}; if isstruct(x.(n))&&isfield(cfg,n), cfg.(n)=merge(cfg.(n),x.(n)); else, cfg.(n)=x.(n); end
end
end
