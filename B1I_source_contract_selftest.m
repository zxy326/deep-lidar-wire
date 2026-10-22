function B1I_source_contract_selftest()
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
cfg=B1I_deep_config(); cfg.deep.dt=.0075;
cfg.deep.initialAttitudeDeg=[25,12,-7]; cfg.deep.armBodyMeters=[1,.2,-.4];
imu=[(23:.0025:23.05)',zeros(21,6)];
a=B1IScheduleProbe(cfg); b=B1IScheduleProbe(cfg);
BDS.B1I.advance_receiver(a,a.endTime,imu);
for t=23.001:.001:23.037, BDS.B1I.advance_receiver(b,t,imu); end
BDS.B1I.advance_receiver(b,b.endTime,imu);
assert(a.ns==5 && isequal(a.updates,b.updates) && isequal(a.queued,b.queued));
assert(a.state==0 && b.state==0 && max(abs(a.updates-(23+(1:5)*.0075)))<1e-12);
short=B1IScheduleProbe(cfg);
BDS.B1I.advance_receiver(short,short.endTime,imu(1:3,:));
assert(short.state==0 && short.ns==0 && ~strcmp(short.stopReason,'Requested interval completed'));
% Code/clock sign and discriminator separation, independently using GPS signs.
c=299792458; chip=c/2.046e6; lambda=c/1561.098e6;
raw=[zeros(2,6),[2e7;2.1e7],[20;-10],ones(2,2)]; sv=raw;
disc=[.02;-.03]; sv(:,7)=raw(:,7)-disc*chip;
rho=raw(:,7)-[5;-6]; rate=[18;-8]; dtr=2e-9; dtv=-1e-10;
[dr,dv,gate]=BDS.B1I.feedback_innovation(raw,sv,rho,rate,dtr,dtv);
gpsCodeChange=-((rho-raw(:,7))/chip+dtr*2.046e6);
gpsCarrierChange=-((rate-raw(:,8))/lambda+dtv*1561.098e6);
assert(max(abs(dr/chip-gpsCodeChange))<1e-12);
assert(max(abs(dv/lambda-gpsCarrierChange))<1e-12);
assert(max(abs(gate-(dr-disc*chip)))<1e-8);
% Nonzero lever arm and nonzero attitude: original source output conventions.
para=BDS.B1I.deep_parameters(cfg,[30,115,100],[3,2,.1]);
f=B1I_navigation_filter(para); input=[.01,-.02,.03,0,0,-9.8];
f.bias(1:3)=[.001,.002,-.001]; before=[f.pos,f.vel,f.P(:)'];
[p,v,sigma]=BDS.B1I.navigation_snapshot(f,input); Cnb=quat2dcm(f.quat);
assert(norm(p-(f.pos+f.arm*Cnb*f.geogInfo.Cn2g))<1e-12);
assert(norm(v-(f.vel+cross(input(1:3)-f.bias(1:3),f.arm)*Cnb))<1e-12);
expected=sqrt(diag(var_phi2angle(f.P(1:3,1:3),Cnb)))'*180/pi;
assert(norm(sigma-expected)<1e-12 && isequal(before,[f.pos,f.vel,f.P(:)']));
% Same-loop diagnostic reuses the source solver, does not affect navigation.
rp=lla2ecef(f.pos); dirs=[1,0,1;-1,0,1;0,1,1;0,-1,1;1,1,.3;-1,1,.2];
rs=rp+dirs./vecnorm(dirs,2,2)*2.2e7;
obs=[rs,zeros(6,3),vecnorm(rs-rp,2,2)+3,zeros(6,1),ones(6,2)];
value=BDS.B1I.source_satnav(obs,rp);
assert(all(isfinite(value)) && norm(value(4:6)-rp)<1e-6 && abs(value(13)*c-3)<1e-6);
assert(norm(value(7:12))<1e-10);
if exist('ecef2lla','file')==2
    reference=satnavSolveWeighted(obs,rp);
    assert(max(abs(value-reference))<1e-6);
end
vs=[1000,2000,500;-1000,200,-400;700,-1200,300;300,2000,-800;900,200,1000;-500,900,600];
Eu=(rs-rp)./vecnorm(rs-rp,2,2); S=sum(Eu.*vs,2); cm=1+S/c;
vp=[3,-2,.5]; drift=1e-9; obs(:,4:6)=vs; obs(:,8)=(S-Eu*vp')./cm+drift*c;
value=BDS.B1I.source_satnav(obs,rp);
assert(norm(value(10:12)-vp)<1e-6 && abs(value(14)-drift)<1e-12);
assert(all(isnan(BDS.B1I.source_satnav(obs(1:3,:),rp))));
for file={'B1I_S_deep.m',fullfile('+BDS','+B1I','deep_receive.m')}
    assert(contains(fileread(fullfile(root,file{1})),'BDS.B1I.advance_receiver('));
end
fprintf('B1I source contracts passed: shared scheduling, raw/filter feedback signs, antenna and Euler outputs, same-loop solver.\n');
end
