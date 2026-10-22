function B1I_source_quality_selftest()
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
cfg=B1I_deep_config(); [p,v,s]=BDS.B1I.quality_options(cfg.deep,0);
assert(~p && v && ~s && ~cfg.deep.sourceCodeVariance && ...
    ~isfield(cfg.deep,'sourceWeightStartIqSeconds'), ...
    'Default must enable rate weighting only, without delayed activation.');
for mask=0:7
    d=struct('sourceRangeWeight',logical(bitget(mask,1)), ...
        'sourceRateWeight',logical(bitget(mask,2)),'sourceCn0Selection',logical(bitget(mask,3)));
    [a,b,c]=BDS.B1I.quality_options(d);
    assert(isequal([a b c],logical(bitget(mask,1:3))));
    d.sourceQualityModel=true; [a,b,c]=BDS.B1I.quality_options(d); assert(a&&b&&c);
end
d=struct('sourceRateWeight',true,'sourceWeightStartIqSeconds',24);
[a,b,c]=BDS.B1I.quality_options(d,23.99); assert(~a&&~b&&~c);
[a,b,c]=BDS.B1I.quality_options(d,24); assert(~a&&b&&~c);
for cn0=[NaN 32.9 33 36.9 37 52]
    [p,v,ip,iv]=BDS.B1I.source_quality(cn0,18,10,10000,90);
    assert(ip==(isfinite(cn0)&&cn0>=33));
    assert(iv==(isfinite(cn0)&&cn0>=37));
    if isfinite(cn0), assert(p>0&&v>0&&isfinite(p)&&isfinite(v)); end
end
[p1,v1]=BDS.B1I.source_quality(40,18,10,10000,45);
e=5*(1-1/(1+exp(-(45-30)/5)));
assert(abs(p1-(sqrt((9e4/.008)/1e4/10)*.5+e)^2)<1e-12);
assert(abs(v1-(.32*18)^3*.0363/1e4*(1575.42/1561.098)^2)<1e-12);
[p2,~]=BDS.B1I.source_quality(40,18,2,500,45);
assert(p2>p1,'Eight integrated milliseconds should be noisier than ten.');
[p3,~]=BDS.B1I.source_quality(40,18,10,10000,10);
assert(p3>p1,'Low-elevation penalty missing.');
[~,~,ip,iv]=BDS.B1I.source_quality(52,18,0,0,45);
assert(~ip&&~iv,'Stale discriminator must not be used.');
fprintf('Source C/N0 thresholds, wavelength scaling, integration time and elevation weighting passed.\n');
end
