function report=B1I_feedback_trial(dataDir)
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
manifest=jsondecode(fileread(fullfile(dataDir,'manifest.json')));
cfg=merge(B1I_deep_config(),manifest.config); cfg.showPlot=false;
cfg.prns=cfg.prns(:)'; cfg.deep.initialAttitudeDeg=cfg.deep.initialAttitudeDeg(:)';
cfg.deep.armBodyMeters=cfg.deep.armBodyMeters(:)';
out=tempname(dataDir); mkdir(out);
diary(fullfile(out,'feedback_ab.log')); cleanup=onCleanup(@()diary('off'));
bcfg=cfg; bcfg.deep.feedbackGain=0;
B=BDS.B1I.deep_receive(bcfg); save(fullfile(out,'control.mat'),'B','-v7.3');
C=BDS.B1I.deep_receive(cfg);
save(fullfile(out,'feedback_ab_inputs.mat'),'B','C','dataDir','-v7.3');
assert(abs(B.initialTimeSeconds-C.initialTimeSeconds)<1e-8);
assert(max(abs(B.diagnostics.epochStartState(1,:)-C.diagnostics.epochStartState(1,:)))<1e-8);
assert(nnz(B.feedbackApplied & (abs(B.codeCorrectionChips)>1e-12|abs(B.carrierCorrectionHz)>1e-12))==0);
[t,ib,ic]=intersect(round(B.timeSeconds(:)*1e6),round(C.timeSeconds(:)*1e6)); t=t/1e6;
truth=readmatrix(fullfile(dataDir,'truth.csv'),'NumHeaderLines',1);
ref=interp1(truth(:,1),truth(:,2:10),t);
report=struct('commonEpochs',numel(t),'firstTime',t(1),'lastTime',t(end), ...
 'note','Same IQ/IMU and initial state. Both runs executed; only feedbackGain changes, including scalar/vector mode. Metrics on common epochs; inspect completion separately.', ...
 'outputDirectory',out);
fig=figure('Position',[100 100 1100 850]); runs={B,C}; ids={ib,ic}; labels={'B no feedback','C feedback'};
for j=1:2
 r=runs{j}; ix=ids{j}; p=zeros(numel(ix),3); v=p;
 for k=1:numel(ix)
  lla=r.positionLla(ix(k),:); p(k,:)=lla2ecef(lla); v(k,:)=r.velocityNed(ix(k),:)*dcmecef2ned(lla(1),lla(2));
 end
 pe=vecnorm(p-ref(:,1:3),2,2); ve=vecnorm(v-ref(:,4:6),2,2); he=r.positionLla(ix,3)-ref(:,9);
 m=struct('positionRmseMeters',sqrt(mean(pe.^2)),'heightRmseMeters',sqrt(mean(he.^2)), ...
  'velocityRmseMps',sqrt(mean(ve.^2)),'maxPositionMeters',max(pe),'maxHeightMeters',max(abs(he)), ...
  'finalHeightMeters',he(end),'updatedEpochs',nnz(r.measurementUsed),'nativeEpochs',numel(r.timeSeconds), ...
  'unlockedChannelEpochs',nnz(r.diagnostics.reason==4),'minimumSatellites',min(r.measurementCount), ...
  'completed',r.timeSeconds(end)>=cfg.startSeconds+cfg.processSeconds-cfg.deep.dt/2,'lastTime',r.timeSeconds(end),'stopReason',r.stopReason,'runDirectory',r.runDir);
 if j==1, report.B=m; else, report.C=m; end
 subplot(3,1,1); plot(t,pe); hold on; ylabel('3D position error (m)'); grid on;
 subplot(3,1,2); plot(t,he); hold on; ylabel('Height error (m)'); grid on;
 subplot(3,1,3); plot(t,ve); hold on; ylabel('3D velocity error (m/s)'); grid on; xlabel('IQ time (s)');
end
for k=1:3, subplot(3,1,k); legend(labels); end
exportgraphics(fig,fullfile(out,'feedback_ab.png')); savefig(fig,fullfile(out,'feedback_ab.fig'));
fid=fopen(fullfile(out,'feedback_ab.json'),'w'); assert(fid>=0); fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true)); fclose(fid);
disp(report.B); disp(report.C); fprintf('OUTPUT: %s\n',out);
end
function cfg=merge(cfg,x)
for name=fieldnames(x)'
 n=name{1}; if isstruct(x.(n))&&isfield(cfg,n), cfg.(n)=merge(cfg.(n),x.(n)); else, cfg.(n)=x.(n); end
end
end

