function report = B1I_simulation_check(dataDirectory,runControl)
% End-to-end B1I IQ / IMU simulation test using the production BB1I_S class.
% Usage: report = B1I_simulation_check;
% Requires MATLAB and a Python with NumPy for input generation (auto-detected).
% Default fixture: 36 seconds, 6 B1I D1 satellites, ~1.18 GB IQ.
% Optional second argument true also runs feedbackGain=0 as an A/B control.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
if nargin<1 || isempty(dataDirectory)
    parent=fullfile(root,'~temp','b1i_simulation');
    if ~isfolder(parent), mkdir(parent); end
    dataDirectory=tempname(parent);
end
if nargin<2, runControl=false; end
validateattributes(runControl,{'logical'},{'scalar'});
if ~isfolder(dataDirectory), mkdir(dataDirectory); end
manifestFile=fullfile(dataDirectory,'manifest.json');
stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
logFile=fullfile(dataDirectory,['validation_',stamp,'.log']);
diary(logFile); logCleanup=onCleanup(@() diary('off'));
report=struct('passed',false,'matlabVersion',version,'dataDirectory',dataDirectory, ...
    'sourceEntry','BDS.B1I.deep_receive -> BB1I_S (same class as B1I_S_deep)', ...
    'errorIdentifier','','errorMessage','');
try
    % Catch component/runtime errors before generating a gigabyte of data.
    fprintf('1/4: MATLAB component self-tests...\n');
    B1I_deep_selftest();
    B1I_receiver_selftest();
    fprintf('2/4: Synchronized signal generation...\n');
    if ~isfile(manifestFile)
        python=find_python();
        generator=fullfile(root,'B1I_simulate.py');
        command=sprintf('"%s" "%s" --output "%s"',python,generator,dataDirectory);
        fprintf('Generating 36 seconds of fc32 IQ (~1.18 GB). Please wait.\n');
        [status,message]=system(command,'-echo');
        assert(status==0,'B1I:SimulationGenerator','Python generator failed: %s',message);
    end
    manifest=jsondecode(fileread(manifestFile));
    cfg=merge(B1I_deep_config(),manifest.config);
    % jsondecode may return column vectors; receiver uses row-vector settings.
    cfg.prns=cfg.prns(:)';
    cfg.deep.initialAttitudeDeg=cfg.deep.initialAttitudeDeg(:)';
    cfg.deep.armBodyMeters=cfg.deep.armBodyMeters(:)';
    cfg.showPlot=false;
    iqInfo=dir(cfg.iqFile);
    assert(isscalar(iqInfo)&&iqInfo.bytes==manifest.iqBytes,'B1I:SimulationSize','IQ fixture length mismatch.');
    truth=readmatrix(fullfile(dataDirectory,'truth.csv'),'NumHeaderLines',1);
    assert(size(truth,2)==13 && all(isfinite(truth(:))),'B1I:Truth','Invalid truth file.');
    % No truth positions, ephemerides or GNSS measurements are supplied to receiver.
    fprintf('3/4: Running B1I acquisition -> tracking -> D1 -> INS -> tracking feedback...\n');
    result=BDS.B1I.deep_receive(cfg);
    fprintf('4/4: Comparing receiver output with independent trajectory truth...\n');
    report=evaluate(result,truth,report,cfg);
    report.receiverResultDirectory=result.runDir;
    report.control=[];
    if runControl
        fprintf('Running A/B control with feedbackGain=0 (same IQ and IMU)...\n');
        cfg.deep.feedbackGain=0;
        control=BDS.B1I.deep_receive(cfg);
        report.control=evaluate(control,truth,struct(),cfg);
        report.control.receiverResultDirectory=control.runDir;
    end
    save(fullfile(dataDirectory,'simulation_validation.mat'),'report','manifest','-v7.3');
    write_report(fullfile(dataDirectory,'simulation_validation.json'),report);
    plot_errors(result,truth,dataDirectory);
    B1I_plot_results(result,truth,fullfile(dataDirectory,'navigation_plots'));
    B1I_height_report(dataDirectory);
    fprintf('\nB1I simulation passed: %d\nPosition RMSE: %.3f m\nVelocity RMSE: %.3f m/s\n', ...
        report.passed,report.positionRmseMeters,report.velocityRmseMps);
    fprintf('Updated epochs: %d/%d; feedback actions: %d; stop: %s\n', ...
        report.updatedEpochs,report.epochs,report.feedbackActions,report.stopReason);
    fprintf('Results and log: %s\n',dataDirectory);
    assert(report.passed,'B1I:SimulationFailed', ...
        'Simulation failed acceptance checks. See simulation_validation.json and validation log.');
catch err
    report.passed=false; report.errorIdentifier=err.identifier; report.errorMessage=err.message;
    write_report(fullfile(dataDirectory,'simulation_validation.json'),report);
    fprintf('\nSimulation stopped: %s\n%s\nLog: %s\n',err.identifier,err.message,logFile);
    rethrow(err);
end
end

function python=find_python()
candidates={fullfile(getenv('USERPROFILE'),'.cache','codex-runtimes', ...
    'codex-primary-runtime','dependencies','python','python.exe')};
try
    env=pyenv;
    if strlength(env.Executable)>0, candidates{end+1}=char(env.Executable); end
catch
    % pyenv is optional; try the command on PATH below.
end
candidates{end+1}='python';
for k=1:numel(candidates)
    [status,~]=system(sprintf('"%s" -c "import numpy"',candidates{k}));
    if status==0, python=candidates{k}; return; end
end
error('B1I:Python','No Python with NumPy found. Configure pyenv or install NumPy in your Python.');
end

function cfg=merge(cfg,override)
names=fieldnames(override);
for j=1:numel(names)
    name=names{j};
    if isstruct(override.(name)) && isfield(cfg,name)
        cfg.(name)=merge(cfg.(name),override.(name));
    else
        cfg.(name)=override.(name);
    end
end
end

function report=evaluate(result,truth,report,cfg)
t=result.timeSeconds; assert(~isempty(t),'B1I:EmptySimulation','No fusion epochs were produced.');
expected=interp1(truth(:,1),truth(:,2:7),t,'linear');
position=zeros(numel(t),3); velocity=position;
for k=1:numel(t)
    position(k,:)=lla2ecef(result.positionLla(k,:));
    Cen=dcmecef2ned(result.positionLla(k,1),result.positionLla(k,2));
    velocity(k,:)=result.velocityNed(k,:)*Cen;
end
dp=sqrt(sum((position-expected(:,1:3)).^2,2));
dv=sqrt(sum((velocity-expected(:,4:6)).^2,2));
truthLla=interp1(truth(:,1),truth(:,8:10),t,'linear');
heightError=result.positionLla(:,3)-truthLla(:,3);
report.heightRmseMeters=sqrt(mean(heightError.^2));
report.maxHeightErrorMeters=max(abs(heightError));
report.initialHeightErrorMeters=heightError(1);
report.finalHeightErrorMeters=heightError(end);
truthDown=interp1(truth(:,1),truth(:,13),t,'linear');
downVelocityError=result.velocityNed(:,3)-truthDown;
report.downVelocityRmseMps=sqrt(mean(downVelocityError.^2));
report.meanHeightErrorMeters=mean(heightError);
report.meanDownVelocityErrorMps=mean(downVelocityError);
report.epochs=numel(t); report.updatedEpochs=nnz(result.measurementUsed);
report.feedbackActions=nnz(result.feedbackApplied & ...
    (abs(result.codeCorrectionChips)>1e-12|abs(result.carrierCorrectionHz)>1e-12));
report.positionRmseMeters=sqrt(mean(dp.^2)); report.maxPositionErrorMeters=max(dp);
report.velocityRmseMps=sqrt(mean(dv.^2)); report.maxVelocityErrorMps=max(dv);
report.firstTime=t(1); report.lastTime=t(end); report.stopReason=result.stopReason;
moving=t>=25;
report.movingEpochs=nnz(moving);
report.movingPositionRmseMeters=sqrt(mean(dp(moving).^2));
report.movingVelocityRmseMps=sqrt(mean(dv(moving).^2));
% Initial engineering acceptance limits, not a claim of research-grade accuracy.
report.criteria=struct('minUpdateFraction',.8,'maxPositionRmseMeters',50, ...
    'maxPositionErrorMeters',100,'maxVelocityRmseMps',1.5,'minMovingSeconds',8);
% Synthetic-fixture engineering regression limits, NOT real-data promises.
report.criteria.maxHeightRmseMeters=.5;
report.criteria.maxHeightErrorMeters=2;
report.criteria.maxDownVelocityRmseMps=.05;
report.passed=all(isfinite([dp;dv])) && result.fusionPerformed && ...
    nnz(result.measurementUsed)/numel(t)>=.8 && t(end)>=cfg.processSeconds-.05 && ...
    nnz(moving)*cfg.deep.dt>=8 && report.positionRmseMeters<50 && ...
report.maxPositionErrorMeters<100 && report.velocityRmseMps<1.5 && ...
    report.heightRmseMeters<report.criteria.maxHeightRmseMeters && ...
    report.maxHeightErrorMeters<report.criteria.maxHeightErrorMeters && ...
    report.downVelocityRmseMps<report.criteria.maxDownVelocityRmseMps && ...
    (cfg.deep.feedbackGain==0 || report.feedbackActions>0);
end

function plot_errors(result,truth,folder)
t=result.timeSeconds; expected=interp1(truth(:,1),truth(:,2:7),t);
position=zeros(numel(t),3); velocity=position;
for k=1:numel(t)
    position(k,:)=lla2ecef(result.positionLla(k,:));
    velocity(k,:)=result.velocityNed(k,:)*dcmecef2ned(result.positionLla(k,1),result.positionLla(k,2));
end
f=figure('Name','B1I simulated truth comparison');
subplot(3,1,1); plot(t,sqrt(sum((position-expected(:,1:3)).^2,2))); ylabel('Position error (m)'); grid on;
subplot(3,1,2); plot(t,sqrt(sum((velocity-expected(:,4:6)).^2,2))); ylabel('Velocity error (m/s)'); grid on;
subplot(3,1,3); stairs(t,result.measurementCount); ylabel('Accepted satellites'); xlabel('IQ file seconds'); grid on;
exportgraphics(f,fullfile(folder,'simulation_errors.png'));
end

function write_report(filename,report)
fid=fopen(filename,'w','n','UTF-8'); assert(fid>=0,'B1I:Report','Cannot write report.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
end
