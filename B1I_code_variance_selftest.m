function B1I_code_variance_selftest
% Monte Carlo check of the E/L covariance formula; no receiver truth aiding.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
previous=rng; cleanup=onCleanup(@()rng(previous)); rng(7411);
cfg=B1I_config; cfg.sampleFreq=4092000; cfg.codeReplica='bandlimited';
model=BDS.B1I.replica(BDS.B1I.codeGene(19),cfg);
chips=(0:4091)/2+.23;
p=BDS.B1I.sample_replica(model,chips); e=BDS.B1I.sample_replica(model,chips+.5); l=BDS.B1I.sample_replica(model,chips-.5);
d=zeros(600,1); v=d;
for k=1:numel(d)
 base=p+2*(randn(size(p))+1i*randn(size(p)));
 E=abs(sum(base.*e)); L=abs(sum(base.*l)); d(k)=.5*(E-L)/(E+L)/model.gain;
 v(k)=BDS.B1I.code_noise_variance(base,p,e,l,model.gain);
end
ratio=var(d)/mean(v);
assert(all(isfinite(v))&&all(v>0)&&ratio>.75&&ratio<1.25,'B1I:CodeVariance','Variance calibration failed: %.4f',ratio);
fprintf('Code variance Monte Carlo passed; measured/predicted = %.4f\n',ratio);
end
