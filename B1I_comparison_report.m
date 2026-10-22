function report = B1I_comparison_report(folder)
% Rebuild comparison plots/metrics from saved runs, without receiver replay.
s=load(fullfile(folder,'comparison_inputs.mat'),'A','B','C','dataDirectory');
truth=readmatrix(fullfile(s.dataDirectory,'truth.csv'),'NumHeaderLines',1);
A=s.A; B=s.B; C=s.C;
assert(isfield(A,'position')&&~isempty(B.timeSeconds)&&~isempty(C.timeSeconds), ...
    'B1I:CompareEmpty','All three runs must provide navigation results.');
assert(abs(B.initialTimeSeconds-C.initialTimeSeconds)<1e-8, ...
    'B1I:CompareInitialization','B and C have different initialization epochs.');
assert(max(abs(B.diagnostics.epochStartState(1,:)-C.diagnostics.epochStartState(1,:)))<1e-8, ...
    'B1I:CompareInitialization','B and C initial navigation states differ.');
[va,validV]=BDS.B1I.scalar_velocity(A);
a=struct('t',A.position.relativeSeconds,'p',A.position.ecef,'v',va, ...
    'validP',A.position.valid,'validV',validV);
b=navigation(B); c=navigation(C);
start=max(B.timeSeconds(1),C.timeSeconds(1));
finish=C.config.startSeconds+C.config.processSeconds;
evalTimes=A.position.relativeSeconds;
evalTimes=evalTimes(evalTimes>=start-1e-8 & evalTimes<=finish+1e-8);
assert(~isempty(evalTimes),'B1I:CompareTime','No shared scalar evaluation epochs.');
[ia,pa]=match(evalTimes,a.t); [ib,pb]=match(evalTimes,b.t); [ic,pc]=match(evalTimes,c.t);
okA=pa&a.validP(ia); okB=pb&b.validP(ib); okC=pc&c.validP(ic);
commonP=okA&okB&okC;
commonV=commonP&a.validV(ia)&b.validV(ib)&c.validV(ic);
report=struct('windowSeconds',[start finish],'expectedScalarEpochs',numel(evalTimes), ...
    'commonPositionEpochs',nnz(commonP),'commonVelocityEpochs',nnz(commonV), ...
    'commonPositionTimes',evalTimes(commonP),'commonVelocityTimes',evalTimes(commonV), ...
    'evaluationNote','Exact common output times, tolerance 1 microsecond. Only truth is interpolated. No interpolated receiver fixes.', ...
    'initializationNote','B/C identical initial navigation states checked; A uses independent scalar positioning with the same IQ/configuration.', ...
    'scope','Single synthetic scenario, not a real-data accuracy certification.', ...
    'directories',struct('A',A.runDir,'B',B.runDir,'C',C.runDir));
report.A=score(a,truth,ia,evalTimes,commonP,commonV,okA,start,finish,A.config.navStepSeconds);
report.B=score(b,truth,ib,evalTimes,commonP,commonV,okB,start,finish,B.config.deep.dt);
report.C=score(c,truth,ic,evalTimes,commonP,commonV,okC,start,finish,C.config.deep.dt);
report.B.stopReason=B.stopReason; report.C.stopReason=C.stopReason;
report.B.completed=B.timeSeconds(end)>=finish-B.config.deep.dt/2;
report.C.completed=C.timeSeconds(end)>=finish-C.config.deep.dt/2;
report.B.updateFraction=mean(B.measurementUsed); report.C.updateFraction=mean(C.measurementUsed);
report.B.nonzeroFeedbackActions=nnz(B.feedbackApplied & (abs(B.codeCorrectionChips)>1e-12|abs(B.carrierCorrectionHz)>1e-12));
report.C.nonzeroFeedbackActions=nnz(C.feedbackApplied & (abs(C.codeCorrectionChips)>1e-12|abs(C.carrierCorrectionHz)>1e-12));
assert(report.B.nonzeroFeedbackActions==0,'B1I:ControlFeedback','Control unexpectedly applied nonzero feedback.');
report.B.inertial=diagnostics(B,s.dataDirectory,folder,'B');
report.C.inertial=diagnostics(C,s.dataDirectory,folder,'C');
report.positionComparisonSufficient=nnz(commonP)>=3;
report.velocityComparisonSufficient=nnz(commonV)>=3;
report.improvement=struct('BminusA_positionRmseMeters',report.B.position3D.rmse-report.A.position3D.rmse, ...
    'CminusB_positionRmseMeters',report.C.position3D.rmse-report.B.position3D.rmse, ...
    'BminusA_velocityRmseMps',report.B.velocity3D.rmse-report.A.velocity3D.rmse, ...
    'CminusB_velocityRmseMps',report.C.velocity3D.rmse-report.B.velocity3D.rmse, ...
    'interpretation','Negative means lower error on common epochs, not statistical significance or improved availability.');
report.metricsNote='RMSE, nearest-rank P95 and maximum over listed common epochs; start-up included. Native availability includes missing tail epochs.';
f=figure('Name','B1I three-way comparison','Position',[100 100 1000 750]);
runs={a,b,c}; ids={ia,ib,ic}; labels={'A scalar','B fusion, no feedback','C fusion + feedback'};
for j=1:3
    expected=interp1(truth(:,1),truth(:,2:7),evalTimes);
    dp=vecnorm(runs{j}.p(ids{j},:)-expected(:,1:3),2,2); dp(~commonP)=NaN;
    dv=vecnorm(runs{j}.v(ids{j},:)-expected(:,4:6),2,2); dv(~commonV)=NaN;
    subplot(3,1,1); plot(evalTimes,dp,'.-'); hold on; ylabel('3D position error (m)'); grid on;
    subplot(3,1,2); plot(evalTimes,dv,'.-'); hold on; ylabel('3D velocity error (m/s)'); grid on;
end
subplot(3,1,1); legend(labels,'Location','best');
subplot(3,1,3); plot(evalTimes,[okA okB okC],'o-'); ylim([-.1 1.1]);
xlim([evalTimes(1)-.1,evalTimes(end)+.1]);
ylabel('Available position (0/1)'); xlabel('IQ seconds'); grid on;
legend(labels,'Location','best');
exportgraphics(f,fullfile(folder,'three_way_errors.png')); savefig(f,fullfile(folder,'three_way_errors.fig'));
save(fullfile(folder,'comparison_report.mat'),'report');
fid=fopen(fullfile(folder,'comparison_report.json'),'w'); assert(fid>=0); cl=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true)); clear cl;
fid=fopen(fullfile(folder,'comparison_summary.txt'),'w'); assert(fid>=0); cl=onCleanup(@() fclose(fid));
fprintf(fid,'Same-data B1I comparison. No accuracy thresholds changed.\n');
fprintf(fid,'Common position epochs: %d/%d; common velocity epochs: %d\n',nnz(commonP),numel(evalTimes),nnz(commonV));
fprintf(fid,'Group | position RMSE (m) | velocity RMSE (m/s) | native position availability\n');
for n={'A','B','C'}
    r=report.(n{1});
    fprintf(fid,'%s | %.6f | %.6f | %.4f\n',n{1},r.position3D.rmse,r.velocity3D.rmse,r.nativePositionAvailability);
end
fprintf(fid,'A has no IMU attitude/bias estimates. B/C attitude and bias metrics are separate.\n');
fprintf(fid,'Small common sample count and one simple trajectory do not establish general superiority.\n');
fprintf(fid,'Inspect JSON for horizontal/vertical, axis errors, P95, maximum and missing-result intervals.\n');
clear cl;
fprintf('%s\n',fileread(fullfile(folder,'comparison_summary.txt')));
fprintf('Comparison report saved: %s\n',folder);
end
function n=navigation(r)
n=struct('t',r.timeSeconds,'p',NaN(numel(r.timeSeconds),3),'v',NaN(numel(r.timeSeconds),3));
for k=1:numel(n.t)
    n.p(k,:)=lla2ecef(r.positionLla(k,:));
    n.v(k,:)=r.velocityNed(k,:)*dcmecef2ned(r.positionLla(k,1),r.positionLla(k,2));
end
n.validP=all(isfinite(n.p),2); n.validV=all(isfinite(n.v),2);
end
function [index,valid]=match(grid,t)
index=ones(size(grid)); valid=false(size(grid));
for k=1:numel(grid)
    [d,index(k)]=min(abs(t-grid(k))); valid(k)=d<1e-6;
end
end
function r=score(n,truth,index,grid,cp,cv,available,start,finish,step)
expected=interp1(truth(:,1),truth(:,2:13),grid);
dp=n.p(index,:)-expected(:,1:3); dv=n.v(index,:)-expected(:,4:6);
pn=dp; vn=dv;
for k=1:numel(grid)
    R=dcmecef2ned(expected(k,7),expected(k,8)); pn(k,:)=dp(k,:)*R'; vn(k,:)=dv(k,:)*R';
end
r=struct('position3D',metric(vecnorm(dp(cp,:),2,2)), ...
    'positionHorizontal',metric(vecnorm(pn(cp,1:2),2,2)), ...
    'positionVertical',metric(abs(pn(cp,3))), ...
    'velocity3D',metric(vecnorm(dv(cv,:),2,2)), ...
    'positionAxisRmseNed',axisrmse(pn(cp,:)), 'velocityAxisRmseNed',axisrmse(vn(cv,:)), ...
    'commonGridAvailability',mean(available));
first=n.t(1)+ceil((start-n.t(1)-1e-8)/step)*step;
native=(first:step:finish)';
[ni,exists]=match(native,n.t); good=exists&n.validP(ni);
r.nativePositionAvailability=mean(good); r.nativeExpectedEpochs=numel(native);
r.nativeVelocityAvailability=mean(exists&n.validV(ni));
r.longestMissingPositionSeconds=longest(~good)*step;
r.lastResultTime=n.t(end);
end
function m=metric(x)
x=x(isfinite(x)); m=struct('count',numel(x),'rmse',NaN,'p95',NaN,'maximum',NaN);
if isempty(x), return; end
x=sort(abs(x)); m.rmse=sqrt(mean(x.^2)); m.p95=x(max(1,ceil(.95*numel(x)))); m.maximum=max(x);
end
function v=axisrmse(x)
if isempty(x), v=[NaN NaN NaN]; else, v=sqrt(mean(x.^2,1)); end
end
function n=longest(flags)
d=diff([false;flags(:);false]); n=max([0;find(d==-1)-find(d==1)]);
end
function summary=diagnostics(r,dataFolder,folder,label)
t=r.timeSeconds; [truth,source]=B1I_inertial_truth(dataFolder,t);
att=mod(r.attitudeDeg-truth(:,1:3)+180,360)-180;
bias=r.bias-truth(:,4:9);
summary=struct('truthSource',source,'attitudeRmseDeg',axisrmse(att), ...
    'attitudeMaxAbsDeg',max(abs(att),[],1),'gyroBiasRmseDegPerSec',axisrmse(bias(:,1:3))*180/pi, ...
    'accelBiasRmseMps2',axisrmse(bias(:,4:6)));
d=r.diagnostics;
assert(size(d.reason,1)==numel(t),'B1I:DiagnosticsLength','Untrimmed diagnostics.');
reasonCount=numel(r.diagnosticLegend.reason);
summary.reasonCounts=zeros(numel(r.prns),reasonCount);
for j=1:numel(r.prns)
    for code=0:reasonCount-1, summary.reasonCounts(j,code+1)=nnz(d.reason(:,j)==code); end
end
summary.prns=r.prns; summary.legend=r.diagnosticLegend;
summary.longestWithoutGnssUpdateSeconds=longest(~r.measurementUsed)*r.config.deep.dt;
f=figure('Name',[label,' attitude and bias error'],'Position',[100 100 1050 800]);
subplot(3,1,1); plot(t,att); yline(0,'k--'); ylabel('YPR error (deg)'); legend('Yaw','Pitch','Roll'); grid on;
subplot(3,1,2); plot(t,bias(:,1:3)*180/pi); yline(0,'k--'); ylabel('Gyro bias error (deg/s)'); legend('X','Y','Z'); grid on;
subplot(3,1,3); plot(t,bias(:,4:6)); yline(0,'k--'); ylabel('Accel bias error (m/s^2)'); xlabel('IQ seconds'); grid on;
exportgraphics(f,fullfile(folder,[label,'_inertial_errors.png'])); savefig(f,fullfile(folder,[label,'_inertial_errors.fig']));
f=figure('Name',[label,' update diagnostics'],'Position',[100 100 1050 800]);
subplot(3,1,1); plot(t,d.rangeGateResidualMeters); ylabel('Range gate residual (m)'); grid on;
subplot(3,1,2); plot(t,d.rateGateResidualMps); ylabel('Rate gate residual (m/s)'); grid on;
change=mod(d.afterCorrection(:,7:9)-d.beforeCorrection(:,7:9)+180,360)-180;
subplot(3,1,3); plot(t,change); ylabel('Euler correction (deg)'); xlabel('IQ seconds'); legend('Yaw','Pitch','Roll'); grid on;
exportgraphics(f,fullfile(folder,[label,'_diagnostics.png'])); savefig(f,fullfile(folder,[label,'_diagnostics.fig']));
% Long-form channel diagnostics, machine-readable; full states remain in MAT.
N=numel(t); C=numel(r.prns);
tableOut=table(repmat(t,C,1),repelem(r.prns(:),N),d.reason(:),d.accepted(:),d.used(:), ...
    d.rangeGateResidualMeters(:),d.rateGateResidualMps(:),d.postRangeResidualMeters(:),d.postRateResidualMps(:), ...
    r.codeCorrectionChips(:),r.carrierCorrectionHz(:), ...
    'VariableNames',{'iq_seconds','prn','reason','accepted','used','range_gate_m','rate_gate_mps','post_range_m','post_rate_mps','feedback_chips','feedback_hz'});
writetable(tableOut,fullfile(folder,[label,'_channel_diagnostics.csv']));
stateTable=table(t,d.epochStartState,d.beforeCorrection,d.afterCorrection,d.biasAfter, ...
    d.filterCorrection,d.updateReason,r.stateSigma,truth,att,bias, ...
    'VariableNames',{'iq_seconds','epoch_start_LLA_VNED_YPR_bias','pre_correction_LLA_VNED_YPR', ...
    'post_correction_LLA_VNED_YPR','bias_after','error_state_correction','update_reason', ...
    'state_sigma','truth_YPR_gyroBias_accelBias','attitude_error_deg','bias_error_rad_s_m_s2'});
writetable(stateTable,fullfile(folder,[label,'_state_diagnostics.csv']));
end
