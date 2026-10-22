function report=B1I_startup_attitude_trial()
% Fixed-configuration paired diagnostic; no receiver tuning or truth injection.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
base=fullfile(root,'~temp','dynamic_week_20261008');
moving=fullfile(root,'~temp','attitude_pair_20261008');
assert(isfile(fullfile(moving,'manifest.json')),'Generate paired attitude fixture first; see B1I_ATTITUDE_DIAGNOSTIC.md.');
baseline=fullfile(base,'matlab_runs','tp0d7c5709_dcea_4ddc_99bb_5c2a15fc4616','deep_result.mat');
assert(isfile(baseline),'Missing frozen baseline result.');
A=BDS.B1I.load_result(baseline);
m=jsondecode(fileread(fullfile(moving,'manifest.json')));
cfg=merge(B1I_deep_config(),m.config); cfg.showPlot=false;
cfg.prns=cfg.prns(:)';cfg.deep.initialAttitudeDeg=cfg.deep.initialAttitudeDeg(:)';
cfg.deep.armBodyMeters=cfg.deep.armBodyMeters(:)';
assert(strcmp(A.config.iqFile,cfg.iqFile) && isequaln(A.config.deep,cfg.deep), ...
    'B1I:PairConfig','Baseline and rotating-body case must have identical IQ and deep settings.');
B=BDS.B1I.deep_receive(cfg);
report=B1I_startup_attitude_report(A,B,base,moving);
end
function cfg=merge(cfg,x)
for name=fieldnames(x)'
    n=name{1}; if isstruct(x.(n))&&isfield(cfg,n), cfg.(n)=merge(cfg.(n),x.(n)); else, cfg.(n)=x.(n); end
end
end
