function B1I_bootstrap_prior_selftest()
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
cfg=B1I_deep_config(); cfg.positionAtmosphere=false; cfg.sampleFreq=4092000;
p=[30,115,100]; rp=lla2ecef(p);
dirs=[1,0,1;-1,0,1;0,1,1;0,-1,1;1,1,.4;-1,1,.3]; dirs=dirs./vecnorm(dirs,2,2);
rs=rp+2.2e7*dirs; sv=[rs,zeros(6,7)]; raw=vecnorm(rs-rp,2,2)+3;
state=rng; clean=onCleanup(@()rng(state)); rng(84721);
boot.channels=cell(1,6);
for j=1:6
 boot.channels{j}=struct('prn',j+18,'boundarySample',(0:2999)'*cfg.sampleFreq/1000, ...
     'codeError',randn(3000,1)*.02,'usable',true(3000,1));
end
[a,clock,C,info]=BDS.B1I.bootstrap_prior(boot,19:24,sv,raw,p+[0,0,.5],2,cfg);
assert(info.accepted && norm(lla2ecef(a)-rp)<1e-4 && abs(clock*299792458-3)<1e-4);
assert(min(eig(C))>0 && max(sqrt(diag(C)))<15);
% Changes in future samples must have no effect on the initialized state/covariance.
for j=1:6, boot.channels{j}.codeError(2001:end)=1e6; end
[a2,c2,C2]=BDS.B1I.bootstrap_prior(boot,19:24,sv,raw,p+[0,0,.5],2,cfg);
assert(isequaln(a,a2)&&isequaln(clock,c2)&&isequaln(C,C2));
% Declared loss of lock must fall back, not manufacture a confident prior.
boot.channels{1}.usable(:)=false;
[~,~,C3,info3]=BDS.B1I.bootstrap_prior(boot,19:24,sv,raw,p,2,cfg);
assert(isempty(C3)&&~info3.accepted);
% Independent scalar loop simulation checks the discrete covariance transfer.
T=.001; wn=cfg.dllBandwidthHz*8*sqrt(.5)/3; a=2*sqrt(.5)*wn; b=wn^2;
phase=0; integrator=0; noise=randn(80000,1); phases=zeros(size(noise)); errors=phases;
for k=1:numel(noise)
    errors(k)=-phase+noise(k);
    integrator=integrator+b*T*errors(k);
    frequency=integrator+a*errors(k);
    phase=phase+T*frequency; phases(k)=phase;
end
A=[1-a*T-b*T*T,-T;b*T,1]; B=[-a*T-b*T*T;b*T]; Q=B*B';
P=reshape((eye(4)-kron(A,A))\Q(:),2,2);
estimated=P(1,1)/(1+P(1,1))*var(errors(10001:end),1);
measured=var(phases(10001:end),1);
assert(measured/estimated>.7 && measured/estimated<1.3,'DLL covariance model failed stochastic loop check.');
fprintf('B1I scalar-prior covariance, causality and unsafe-history fallback tests passed.\n');
end
