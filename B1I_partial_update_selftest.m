function B1I_partial_update_selftest()
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
cfg=B1I_deep_config(); cfg.deep.initialAttitudeDeg=[0,0,0]; cfg.deep.armBodyMeters=[0,0,0];
[ip,iv,reason]=BDS.B1I.observation_masks([1;101;NaN;2;Inf],[1;2;3;16;Inf],cfg.deep);
assert(isequal(ip,[true;false;false;true;false]));
assert(isequal(iv,[true;true;true;false;false]));
assert(isequal(reason,uint8([0;9;9;10;8])));
p=[30,115,100]; rp=lla2ecef(p); rs=rp+[2e7,1e7,1.5e7]; vs=[1000,2000,-500];
for physical=[false,true]
 cfg.deep.physicalRangeRate=physical; para=BDS.B1I.deep_parameters(cfg,p,[0,0,0]);
 for mode=1:3
  f=B1I_navigation_filter(para); coast=B1I_navigation_filter(para);
  for item={f,coast}
   x=item{1}; x.motion.state=1; x.motion.state0=1; x.motion.wN1=Inf; x.motion.vN1=Inf;
  end
  input=[f.geogInfo.wien,0,0,-f.geogInfo.g];
  [rho,rate,u]=rho_rhodot_cal_ecef(rs,vs,rp,[0,0,0]);
  if physical, scale=BDS.B1I.rate_scale(rs,vs,rp); else, scale=1-sum(u.*vs,2)/299792458; end
  obs=[rs,vs,rho+.1,rate/scale+.01,225,1];
  indexP=mode~=2; indexV=mode~=1;
  if ~indexP, obs(7)=NaN; end
  if ~indexV, obs(8)=NaN; end
  f.run(input,obs,indexP,indexV); coast.run(input);
  assert(all(isfinite(f.P(:)))&&min(eig((f.P+f.P')/2))>-1e-8);
  if indexP, assert(trace(f.P(7:10,7:10))<trace(coast.P(7:10,7:10))); end
  if indexV, assert(trace(f.P([4:6,11],[4:6,11]))<trace(coast.P([4:6,11],[4:6,11]))); end
 end
end
fprintf('B1I separate gates and single-satellite P-only/V-only/combined update tests passed.\n');
end
