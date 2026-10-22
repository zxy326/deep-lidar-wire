function report=B1I_source_quality_trial(dataDir)
% Compare the user's saved baseline with source-quality trial on identical IQ/IMU.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
B1I_source_quality_selftest();
s=jsondecode(fileread(fullfile(dataDir,'simulation_validation.json')));
B=BDS.B1I.load_result(fullfile(s.receiverResultDirectory,'deep_result.mat'));
assert(~isfield(B.config.deep,'sourceQualityModel')||~B.config.deep.sourceQualityModel);
cfg=B.config; cfg.showPlot=false; cfg.deep.sourceQualityModel=true;
cfg.deep.sourceCodeVariance=false;
out=tempname(dataDir); mkdir(out);
diary(fullfile(out,'trial.log')); cleanup=onCleanup(@()diary('off'));
C=BDS.B1I.deep_receive(cfg);
assert(abs(B.initialTimeSeconds-C.initialTimeSeconds)<1e-9);
assert(max(abs(B.diagnostics.epochStartState(1,:)-C.diagnostics.epochStartState(1,:)))<1e-9, ...
    'Trial must start from the same navigation/filter state.');
truth=readmatrix(fullfile(dataDir,'truth.csv'),'NumHeaderLines',1);
report=struct('baseline',metrics(B,truth),'sourceQuality',metrics(C,truth), ...
    'outputDirectory',out,'sameInitialState',true, ...
    'note','Saved baseline versus identical IQ/IMU replay; only source C/N0 selection and range/rate noise model changed. No truth enters receiver.');
f=figure('Name','Source quality model A/B','Position',[100 100 1000 850]);
for j=1:2
    if j==1, r=B; else, r=C; end
    ref=interp1(truth(:,1),truth(:,8:13),r.timeSeconds);
    subplot(3,1,1); plot(r.timeSeconds,r.positionLla(:,3)-ref(:,3)); hold on; grid on; ylabel('Height error (m)');
    subplot(3,1,2); plot(r.timeSeconds,1000*(r.velocityNed(:,3)-ref(:,6))); hold on; grid on; ylabel('Down error (mm/s)');
    subplot(3,1,3); stairs(r.timeSeconds,r.measurementCount); hold on; grid on; ylabel('Accepted channels'); xlabel('IQ seconds');
end
for j=1:3, subplot(3,1,j); legend('Current baseline','Source quality trial'); end
exportgraphics(f,fullfile(out,'comparison.png')); savefig(f,fullfile(out,'comparison.fig'));
fid=fopen(fullfile(out,'report.json'),'w'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
disp(report.baseline); disp(report.sourceQuality); fprintf('TRIAL REPORT: %s\n',out);
end
function m=metrics(r,truth)
t=r.timeSeconds; ref=interp1(truth(:,1),truth(:,8:13),t);
h=r.positionLla(:,3)-ref(:,3); d=r.velocityNed(:,3)-ref(:,6);
dp=zeros(numel(t),3);
for k=1:numel(t), dp(k,:)=lla2ecef(r.positionLla(k,:))-lla2ecef(ref(k,1:3)); end
dv=r.velocityNed-ref(:,4:6);
m=struct('runDirectory',r.runDir,'epochs',numel(t),'updates',nnz(r.measurementUsed), ...
    'completed',t(end)>=r.config.startSeconds+r.config.processSeconds-.005, ...
    'stopReason',r.stopReason,'heightRmse',sqrt(mean(h.^2)), ...
    'downRmse',sqrt(mean(d.^2)),'positionRmse',sqrt(mean(sum(dp.^2,2))), ...
    'velocityRmse',sqrt(mean(sum(dv.^2,2))),'finalHeight',h(end));
for item={'startup','stable'}
    name=item{1}; ix=t<=r.initialTimeSeconds+1;
    if strcmp(name,'stable'), ix=t>=r.initialTimeSeconds+2; end
    m.(name)=struct('epochs',nnz(ix),'heightRmse',sqrt(mean(h(ix).^2)), ...
        'heightMean',mean(h(ix)),'heightMaxAbs',max(abs(h(ix))), ...
        'downRmse',sqrt(mean(d(ix).^2)),'downMean',mean(d(ix)), ...
        'positionRmse',sqrt(mean(sum(dp(ix,:).^2,2))), ...
        'velocityRmse',sqrt(mean(sum(dv(ix,:).^2,2))));
end
end
