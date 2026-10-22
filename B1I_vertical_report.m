function stats=B1I_vertical_report(runDirectory,truthFile)
% Truth is used ONLY for post-run validation, never for receiver correction.
r=BDS.B1I.load_result(fullfile(runDirectory,'deep_result.mat'));
truth=readmatrix(truthFile,'NumHeaderLines',1); t=r.timeSeconds(:);
ref=interp1(truth(:,1),truth(:,8:13),t,'linear',NaN);
dh=r.positionLla(:,3)-ref(:,3); dv=r.velocityNed(:,3)-ref(:,6);
assert(~isempty(t)&&all(isfinite([dh;dv])),'B1I:TruthCoverage','Missing finite vertical truth.');
stats=struct('height',measure(dh),'downVelocity',measure(dv), ...
 'firstTime',t(1),'lastTime',t(end),'epochs',numel(t),'updatedEpochs',nnz(r.measurementUsed));
% Fixed elapsed-time late window, distinct from full-window metrics.
late=t>=t(1)+7;
if any(late), stats.lateHeight=measure(dh(late)); stats.lateDownVelocity=measure(dv(late)); end
f=figure('Name','B1I vertical validation','Position',[100 100 1100 850]);
subplot(3,1,1); plot(t,r.positionLla(:,3),t,ref(:,3),'--'); ylabel('Height (m)'); legend('Estimate','Truth'); grid on;
subplot(3,1,2); plot(t,dh); yline(0,'k--'); ylabel('Height error (m)'); grid on;
subplot(3,1,3); plot(t,dv); yline(0,'k--'); ylabel('Down velocity error (m/s)'); xlabel('IQ seconds'); grid on;
exportgraphics(f,fullfile(runDirectory,'vertical_validation.png')); savefig(f,fullfile(runDirectory,'vertical_validation.fig'));
fid=fopen(fullfile(runDirectory,'vertical_validation.json'),'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(stats,'PrettyPrint',true)); disp(stats);
end
function m=measure(e)
a=sort(abs(e)); m=struct('rmse',sqrt(mean(e.^2)),'mean',mean(e),'standardDeviation',std(e), ...
 'p95Absolute',a(max(1,ceil(.95*numel(a)))),'maxAbsolute',max(a),'first',e(1),'last',e(end));
end
