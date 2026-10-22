function report=B1I_startup_attitude_report(A,B,base,moving)
% Error and model-based 3-sigma diagnostics; not a statistical certification.
out=fullfile(moving,'diagnostics'); if ~isfolder(out),mkdir(out);end
assert(isequal(A.timeSeconds,B.timeSeconds),'Incomplete or unequal paired run intervals.');
report=struct('baseline',score(A,base),'rotating',score(B,moving), ...
    'note','One noise seed, identical IQ/translation/settings, different body attitude; 3-sigma is filter-reported uncertainty, not guaranteed confidence.');
truth=readmatrix(fullfile(moving,'truth.csv'),'NumHeaderLines',1);
[inertial,~]=B1I_inertial_truth(moving,B.timeSeconds);
err=[mod(B.attitudeDeg-inertial(:,1:3)+180,360)-180, ...
    (B.bias(:,1:3)-inertial(:,4:6))*180/pi,B.bias(:,4:6)-inertial(:,7:9)];
sigma=[B.attitudeSigmaDeg,B.stateSigma(:,12:14)*180/pi,B.stateSigma(:,15:17)];
labels={{'Yaw (deg)','Pitch (deg)','Roll (deg)'}, ...
    {'Gyro X (deg/s)','Gyro Y (deg/s)','Gyro Z (deg/s)'}, ...
    {'Accel X (m/s^2)','Accel Y (m/s^2)','Accel Z (m/s^2)'}};
names={'attitude_error','gyro_bias_error','accel_bias_error'};
for g=1:3
    f=figure('Name',names{g});
    for j=1:3
        q=(g-1)*3+j;subplot(3,1,j);
        plot(B.timeSeconds,err(:,q),'b',B.timeSeconds,3*sigma(:,q),'r--',B.timeSeconds,-3*sigma(:,q),'r--');
        ylabel(labels{g}{j});grid on;
        if j==1,legend('Estimate - truth','Filter +3 sigma','Filter -3 sigma');end
    end
    xlabel('IQ seconds');savefig(f,fullfile(out,[names{g},'.fig']));exportgraphics(f,fullfile(out,[names{g},'.png']));
end
f=figure('Name','Startup handover');mask=B.timeSeconds<=B.initialTimeSeconds+2;
t=B.timeSeconds(mask);h=interp1(truth(:,1),truth(:,10),t);
subplot(3,1,1);plot(t,A.positionLla(mask,3)-h,t,B.positionLla(mask,3)-h);ylabel('Height error (m)');legend('Fixed body','Rotating body');grid on;
subplot(3,1,2);plot(t,B.diagnostics.beforeCorrection(mask,3)-h,t,B.diagnostics.afterCorrection(mask,3)-h);ylabel('Height error (m)');legend('Before correction','After correction');grid on;
subplot(3,1,3);plot(t,B.diagnostics.afterCorrection(mask,3)-B.diagnostics.beforeCorrection(mask,3));ylabel('Height update (m)');xlabel('IQ seconds');grid on;
savefig(f,fullfile(out,'startup.fig'));exportgraphics(f,fullfile(out,'startup.png'));
save(fullfile(out,'paired_results.mat'),'A','B','report','-v7.3');
fid=fopen(fullfile(out,'report.json'),'w');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));disp(report.baseline);disp(report.rotating);
fprintf('Diagnostic report: %s\n',out);
end
function s=score(r,folder)
t=r.timeSeconds; truth=readmatrix(fullfile(folder,'truth.csv'),'NumHeaderLines',1);
ref=interp1(truth(:,1),truth(:,8:13),t);[it,~]=B1I_inertial_truth(folder,t);
ae=mod(r.attitudeDeg-it(:,1:3)+180,360)-180;
be=r.bias-it(:,4:9); sg=[r.stateSigma(:,12:14),r.stateSigma(:,15:17)];
startup=t<=r.initialTimeSeconds+2;tail=t>=t(end)-10; early=t<24;
dh=r.positionLla(:,3)-ref(:,3);
[~,peakIndex]=max(abs(dh(startup))); startupTimes=t(startup);
s=struct('runDirectory',r.runDir,'epochs',numel(t),'updates',nnz(r.measurementUsed), ...
    'startupHeightPeakMeters',max(abs(dh(startup))), ...
    'startupPeakIqSeconds',startupTimes(peakIndex), ...
    'preMotionHeightPeakMeters',max(abs(dh(early))), ...
    'largestStartupPropagationHeightStepMeters',max(abs(r.diagnostics.beforeCorrection(startup,3)-r.diagnostics.epochStartState(startup,3))), ...
    'largestStartupHeightUpdateMeters',max(abs(r.diagnostics.afterCorrection(startup,3)-r.diagnostics.beforeCorrection(startup,3))), ...
    'tailHeightMeanMeters',mean(dh(tail)),'heightRmseMeters',sqrt(mean(dh.^2)), ...
    'velocityRmseMps',sqrt(mean(sum((r.velocityNed-ref(:,4:6)).^2,2))), ...
    'attitudeRmseDeg',sqrt(mean(ae.^2)),'tailAttitudeRmseDeg',sqrt(mean(ae(tail,:).^2)), ...
    'finalBiasError',be(end,:),'tailBiasRmse',sqrt(mean(be(tail,:).^2)), ...
    'attitudeWithin3SigmaFraction',mean(abs(ae)<=3*r.attitudeSigmaDeg), ...
    'biasWithin3SigmaFraction',mean(abs(be)<=3*sg));
end
