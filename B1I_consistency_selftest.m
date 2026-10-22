function B1I_consistency_selftest()
% Regression: parameter wiring, legacy files and initialization range rates.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
cfg=B1I_deep_config(); cfg.deep.initialAttitudeDeg=[0,0,0];
cfg.deep.armBodyMeters=[0,0,0]; position=[45,126,100];
for enabled=[false,true]
    cfg.deep.physicalRangeRate=enabled;
    para=BDS.B1I.deep_parameters(cfg,position,[1,2,3]);
    f=B1I_navigation_filter(para);
    assert(f.physicalRangeRate==enabled);
end
% Both public entry points must use the same parameter factory.
for file={'B1I_S_deep.m',fullfile('+BDS','+B1I','deep_receive.m')}
    assert(contains(fileread(fullfile(root,file{1})),'BDS.B1I.deep_parameters('));
end
file=[tempname,'.mat']; cleanup=onCleanup(@()remove_file(file));
result=struct('timeSeconds',[1;2],'positionLla',ones(2,3));
save(file,'result'); assert(isequaln(BDS.B1I.load_result(file),result));
b1iDeepResult=result; save(file,'b1iDeepResult');
assert(isequaln(BDS.B1I.load_result(file),result));
% Deterministic full-rank geometry; independently form received Doppler rates.
rp=lla2ecef(position); Cen=dcmecef2ned(position(1),position(2));
directions=[1,0,1;-1,0,1;0,1,1;0,-1,1;1,1,.3;-1,1,.2];
directions=directions./vecnorm(directions,2,2);
rs=rp+2.2e7*directions; vs=[1000,2200,400;-1800,500,800; ...
    700,-1800,900;-900,-1500,400;2200,600,-400;1300,-800,1100];
v=[3,-2,.4]; vp=v*Cen; u=(rp-rs)./vecnorm(rp-rs,2,2);
clock=2e-9; geometric=sum(u.*(vp-vs),2);
for enabled=[false,true]
    cfg.deep.physicalRangeRate=enabled;
    s=ones(6,1);
    if enabled
        rotationDerivative=[7.292115e-5*rs(:,2),-7.292115e-5*rs(:,1),zeros(6,1)];
        s=1+sum(u.*(rotationDerivative-vs),2)/299792458;
    end
    sv=[rs,vs,zeros(6,1),geometric./s+clock*299792458,ones(6,2)];
    cfg.deep.initialVelocityNed=[NaN,NaN,NaN];
    [estimated,drift]=BDS.B1I.initial_velocity(sv,position,cfg);
    assert(norm(estimated-v)<1e-8 && abs(drift-clock)<1e-15);
    cfg.deep.initialVelocityNed=v;
    [estimated,drift]=BDS.B1I.initial_velocity(sv,position,cfg);
    assert(isequal(estimated,v)&&abs(drift-clock)<1e-15);
end
fprintf('B1I entry/schema/initial-velocity consistency tests passed.\n');
end
function remove_file(file)
if isfile(file), delete(file); end
end
