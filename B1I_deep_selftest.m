function B1I_deep_selftest()
% MATLAB regression tests. Synthetic fixtures only, NOT a real-data claim.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
cfg=B1I_deep_config(); cfg.showPlot=false;
% Fractional-phase replica and calibrated discriminator, independent of INS.
waveCfg=cfg; waveCfg.sampleFreq=4092000; waveCfg.codeReplica='bandlimited';
model=BDS.B1I.replica(BDS.B1I.codeGene(19),waveCfg);
for phase=[.03 .23 .48]
    chips=(0:4091)/2+phase;
    for delta=[-.04 0 .04]
        xWave=BDS.B1I.sample_replica(model,chips+delta);
        E=abs(sum(xWave.*BDS.B1I.sample_replica(model,chips+.5)));
        L=abs(sum(xWave.*BDS.B1I.sample_replica(model,chips-.5)));
        de=.5*(E-L)/(E+L)/model.gain;
        assert(abs(de-delta)<.008,'Fractional code discriminator sign/gain regression.');
    end
end
expect_error(@() BDS.B1I.load_imu(cfg),'B1I:SyncRequired');
% Test explicit time mapping, units and proper sensor -> FRD body rotation.
imuPath=[tempname,'.mat']; iqPath=[tempname,'.dat'];
cleanup=onCleanup(@() clean_files(imuPath,iqPath));
imu=[(100:0.01:101)',ones(101,1)*180,zeros(101,2),zeros(101,2),-ones(101,1)];
save(imuPath,'imu'); cfg.imuFile=imuPath; cfg.imu.format='mat';
cfg.sync.confirmed=true; cfg.sync.imuReferenceSeconds=100;
cfg.sync.imuToIqOffsetSeconds=30; cfg.sync.imuClockScale=1.00001;
cfg.sync.imuAxesToBody=[0,-1,0;1,0,0;0,0,1];
cfg.imu.specificForceConfirmed=true; cfg.imu.gyroUnits='deg/s'; cfg.imu.accelUnits='g';
converted=BDS.B1I.load_imu(cfg);
assert(abs(converted(end,1)-31.00001)<1e-10);
assert(abs(converted(1,3)-pi)<1e-12 && abs(converted(1,7)+9.80665)<1e-12);
bad=cfg; bad.sync.imuAxesToBody=diag([1,1,-1]);
expect_error(@() BDS.B1I.load_imu(bad),'B1I:Axes');
bad=cfg; bad.imu.specificForceConfirmed=false;
expect_error(@() BDS.B1I.load_imu(bad),'B1I:SpecificForce');
% Continuous and segmented tracking must preserve exactly the same state.
cfg.iqFile=iqPath; cfg.startSeconds=0; cfg.processSeconds=1.05;
fs=cfg.sampleFreq; samples=0:round(fs*.75)-1; frequency=1850; delay=1234;
chips=(samples-delay)*2.046e6*(1+frequency/1561.098e6)/fs;
code=BDS.B1I.codeGene(19);
x=code(mod(floor(chips),2046)+1).*(-1).^floor(chips/2046).* ...
    exp(2i*pi*(frequency*samples/fs+.18));
x=[x,zeros(1,round(.3*fs))]; % Actual loss of signal must still stop the channel.
fid=fopen(iqPath,'wb','ieee-le'); assert(fid>=0);
closer=onCleanup(@() fclose(fid)); fwrite(fid,[real(x);imag(x)],'single'); clear closer;
acq=table(19,delay,1750,'VariableNames',{'PRN','CodeDelaySamples','DopplerHz'});
a=BDS.B1I.DeepChannel(cfg,acq); b=BDS.B1I.DeepChannel(cfg,acq);
channelsCleanup=onCleanup(@() close_channels(a,b));
a.advance(.7);
for t=.01:.01:.7, b.advance(t); end
assert(a.n==b.n && a.ptr==b.ptr && abs(a.prompt-b.prompt)<1e-8);
assert(abs(a.carrNco-b.carrNco)<1e-10 && abs(a.remCode-b.remCode)<1e-10);
assert(a.locked && abs(a.carrNco-cfg.ifHz-frequency)<3);
scalarCfg=cfg; scalarCfg.processSeconds=.7;
scalar=BDS.B1I.track(scalarCfg,acq);
assert(a.n==numel(scalar.prompt) && abs(a.prompt-scalar.prompt(end))<1e-8);
% Feedback sign, bound and actual effect on subsequent correlations.
anchor=struct('epoch',1,'sow',100000);
[before,~,ok]=a.observe(.7,anchor); assert(ok);
oldCode=a.remCode; oldCarrier=a.carrNco;
a.feedback(100,10);
[after,~,ok]=a.observe(.7,anchor); assert(ok && after>before);
assert(a.remCode>oldCode && a.carrNco>oldCarrier);
assert(a.remCode-oldCode<=cfg.deep.maxCodeCorrectionChips+1e-12);
assert(a.carrNco-oldCarrier<=cfg.deep.maxCarrierCorrectionHz+1e-12);
a.advance(.71); b.advance(.71);
assert(abs(a.prompt-b.prompt)>1e-6,'NCO feedback did not affect future IQ processing.');
% Vector switching must disable the independent DLL and separate PLL/NCO.
a.beginFusion(); assert(a.vectorActive && a.codeInt==0);
estimate=a.carrInt; a.feedback(10,1);
assert(a.carrInt==estimate,'Navigation feedback corrupted PLL frequency estimate.');
a.advance(.72);
assert(a.codeInt==0 && a.codeDiscCount>0);
assert(abs(a.codeNco-2.046e6*(1+(a.carrInt-cfg.ifHz)/1561.098e6))<1e-8);
a.resetDiscriminators(); assert(a.codeDiscCount==0);
% One-block vectorPLL2 equation must match the source control law exactly.
oldPtr=a.ptr; oldPhase=a.phase; oldEstimate=a.carrInt; oldNco=a.carrNco;
block=ceil((2046-a.remCode)/a.codeNco*fs); T=block/fs;
a.advance((oldPtr+block+.1)/fs);
assert(a.ptr==oldPtr+block);
wn=cfg.pllBandwidthHz*8*sqrt(.5)/3;
df=oldEstimate-oldNco; dp=a.phaseError-df*T;
assert(abs(a.carrInt-(oldEstimate+wn^2*T*dp))<1e-8);
expectedPhase=mod(oldPhase+oldNco*T+df*T+2*sqrt(.5)*wn*T*dp,1);
assert(abs(mod(a.phase-expectedPhase+.5,1)-.5)<1e-8);
expect_error(@() a.feedback(NaN,0),'B1I:Feedback');
% Acceleration aid must advance BOTH oscillator and estimate at sample time.
a.carrierRate=2.5;
oldPtr=a.ptr; oldPhase=a.phase; oldEstimate=a.carrInt; oldNco=a.carrNco;
block=ceil((2046-a.remCode)/a.codeNco*fs); T=block/fs;
a.advance((oldPtr+block+.1)/fs);
dp=a.phaseError-(oldEstimate-oldNco)*T;
assert(abs(a.carrInt-(oldEstimate+2.5*T+wn^2*T*dp))<1e-8);
assert(abs(a.carrNco-(oldNco+2.5*T))<1e-8);
expectedPhase=mod(oldPhase+oldEstimate*T+.5*2.5*T^2+2*sqrt(.5)*wn*T*dp,1);
assert(abs(mod(a.phase-expectedPhase+.5,1)-.5)<1e-8);
[~,d1]=a.observe(a.ptr/fs,anchor);
[~,d2]=a.observe(a.ptr/fs+.0004,anchor);
assert(abs(d2-d1-.001)<1e-9,'Doppler epoch transport failed.');
a.carrierRate=0;
b.advance(1.05);
assert(b.stopped && ~b.locked,'True sustained signal loss must stop the channel.');
% Close owned readers BEFORE attempting Windows fixture deletion.
clear channelsCleanup;
% Original INS/filter interface: stationary, known attitude and gyro Earth rate.
cfg.deep.initialAttitudeDeg=[0,0,0]; cfg.deep.armBodyMeters=[0,0,0];
p=[30,115,100]; para=BDS.B1I.deep_parameters(cfg,p,[0,0,0]);
f=filter_single(para); f.motion.state=1; f.motion.wN1=Inf; f.motion.vN1=Inf;
input=[f.geogInfo.wien,0,0,-f.geogInfo.g];
rp=lla2ecef(p); Cen=dcmecef2ned(p(1),p(2));
az=(0:60:300)'; el=[30;45;60;35;55;70];
losN=[cosd(el).*cosd(az),cosd(el).*sind(az),-sind(el)];
rs=rp+2.3e7*(losN*Cen); rho=sqrt(sum((rs-rp).^2,2));
sv=[rs,zeros(6,3),rho,zeros(6,1),ones(6,1)*225,ones(6,1)];
% Known inertial acceleration projects to Doppler rate with the correct sign.
probe=filter_single(para); probe.wdot=[0,0,0];
accelNed=[.4,.25,0];
probeInput=[probe.geogInfo.wien,accelNed+[0,0,-probe.geogInfo.g]];
predRate=BDS.B1I.carrier_rate(probe,probeInput,sv,zeros(6,4),rp,[0,0,0]);
[~,~,los]=rho_rhodot_cal_ecef(rs,zeros(6,3),rp,[0,0,0]);
expectedRate=-los*(accelNed*Cen)'/(299792458/1561.098e6);
assert(max(abs(predRate-expectedRate))<1e-8,'INS acceleration/Doppler sign regression.');
logged=B1I_diagnostic_filter(para);
legacyPara=para; legacyPara.physicalRangeRate=false;
legacy=B1I_navigation_filter(legacyPara);
legacy.motion.state=1; legacy.motion.wN1=Inf; legacy.motion.vN1=Inf;
logged.motion.state=1; logged.motion.wN1=Inf; logged.motion.vN1=Inf;
for j=1:100
    f.run(input,sv,true(6,1),true(6,1));
    logged.run(input,sv,true(6,1),true(6,1));
    legacy.run(input,sv,true(6,1),true(6,1));
end
assert(isequaln(f.P,logged.P)&&isequaln(f.quat,logged.quat)&& ...
    isequaln(f.pos,logged.pos)&&isequaln(f.vel,logged.vel)&&isequaln(f.bias,logged.bias), ...
    'Diagnostic instrumentation changed filter numerics.');
assert(numel(logged.beforeCorrection)==9&&numel(logged.correctionState)==9);
assert(isequaln(f.P,legacy.P)&&isequaln(f.vel,legacy.vel)&&isequaln(f.bias,legacy.bias), ...
    'Disabling physical range-rate must retain original filter numerics.');
% Independent finite-difference light-time check for a moving satellite.
sat=[2e7,1.5e7,1e7]; satvel=[-1000,2400,700]; recvvel=[2,1,.5];
[d0,rotpos,rotvel]=flight(0,sat,satvel,rp,recvvel);
h=.01; numeric=(flight(h,sat,satvel,rp,recvvel)-flight(-h,sat,satvel,rp,recvvel))/(2*h);
[~,geometric]=rho_rhodot_cal_ecef(rotpos,rotvel,rp,recvvel);
physical=geometric/BDS.B1I.rate_scale(rotpos,rotvel,rp);
assert(d0>0 && abs(physical-numeric)<1e-4,'Light-time range-rate regression.');
assert(norm(lla2ecef(f.pos)-rp)<.1 && norm(f.vel)<.1);
assert(all(isfinite(f.P(:))) && norm(f.P-f.P','fro')<1e-8);
assert(min(eig(f.P))>-1e-8);
% Receiver class contract and checked endpoint interpolation (no real IQ read).
r=BB1I_S(cfg); receiverCleanup=onCleanup(@() delete(r));
assert(r.state==0 && ~r.initialized && isnan(r.tp));
expect_error(@() r.run(.001),'B1I:State');
expect_error(@() r.imu_input(0,zeros(1,6)),'B1I:State');
expect_error(@() r.channel_vector(),'B1I:State');
required={'initialize','run','imu_input','channel_vector','get_result', ...
    'clean_storage','save_ephemeris','plot_result','stop'};
assert(all(ismember(required,methods(r))));
rows=[0,zeros(1,6);.01,ones(1,6);.02,2*ones(1,6)];
assert(max(abs(BDS.B1I.imu_at(rows,.005,cfg)-.5))<1e-12);
expect_error(@() BDS.B1I.imu_at(rows,.03,cfg),'B1I:IMUCoverage');
rows(3,1)=.1;
expect_error(@() BDS.B1I.imu_at(rows,.05,cfg),'B1I:IMUGap');
fprintf('B1I deep synthetic component tests passed; real synchronized replay is still required.\n');
clear receiverCleanup cleanup;
end
function expect_error(fn,id)
try
    fn();
catch err
    assert(strcmp(err.identifier,id),'Unexpected error: %s',err.message); return;
end
error('B1I:Test','Expected error %s was not raised.',id);
end
function close_channels(a,b)
delete(a); delete(b);
end
function clean_files(varargin)
% Channel cleanup is explicitly completed before this callback on success.
for j=1:nargin, if isfile(varargin{j}), delete(varargin{j}); end, end
end
function [rho,rotpos,rotvel]=flight(t,rs,vs,rp,vp)
tau=.08;
for k=1:10
 theta=7.292115e-5*tau; R=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
 rotpos=(rs+vs*(t-tau))*R'; rotvel=vs*R';
 rho=norm(rp+vp*t-rotpos); tau=rho/299792458;
end
end
