function report = B1I_compare_runs(dataDirectory)
% Run A scalar, B fusion without feedback, C fusion with feedback on same IQ.
% Each invocation creates a unique folder. No previous reports overwritten.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
manifest=jsondecode(fileread(fullfile(dataDirectory,'manifest.json')));
cfg=merge(B1I_deep_config(),manifest.config); cfg.showPlot=false;
cfg.prns=cfg.prns(:)'; cfg.deep.initialAttitudeDeg=cfg.deep.initialAttitudeDeg(:)';
cfg.deep.armBodyMeters=cfg.deep.armBodyMeters(:)';
parent=fullfile(dataDirectory,'comparisons'); if ~isfolder(parent), mkdir(parent); end
folder=tempname(parent); mkdir(folder);
diary(fullfile(folder,'comparison.log')); cleanup=onCleanup(@() diary('off'));
fprintf('Comparison output: %s\n',folder);
try
    B1I_deep_selftest(); B1I_receiver_selftest();
    fprintf('A/3: Full-duration scalar B1I (no IMU)...\n');
    A=BDS.B1I.receive(cfg);
    fprintf('B/3: B1I + IMU, feedback gain zero...\n');
    bcfg=cfg; bcfg.deep.feedbackGain=0; B=BDS.B1I.deep_receive(bcfg);
    fprintf('C/3: B1I + IMU, configured feedback...\n');
    C=BDS.B1I.deep_receive(cfg);
    save(fullfile(folder,'comparison_inputs.mat'),'A','B','C','dataDirectory','-v7.3');
    report=B1I_comparison_report(folder);
catch err
    fid=fopen(fullfile(folder,'failure.txt'),'w');
    if fid>=0, fprintf(fid,'%s\n%s',err.identifier,getReport(err,'extended','hyperlinks','off')); fclose(fid); end
    rethrow(err);
end
end
function cfg=merge(cfg,override)
names=fieldnames(override);
for k=1:numel(names)
    n=names{k};
    if isstruct(override.(n))&&isfield(cfg,n), cfg.(n)=merge(cfg.(n),override.(n));
    else, cfg.(n)=override.(n); end
end
end
