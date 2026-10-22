function figures = B1I_plot_results(source,truth,outputDirectory)
% GPS-style navigation diagnostics from B1I saved results; no receiver rerun.
% B1I_plot_results(simulationDirectory) resolves simulation_validation.json.
% B1I_plot_results(result,truth,outputDirectory) accepts the receiver struct.
% Truth columns: IQ seconds, ECEF position/velocity, LLA, NED velocity.
if nargin<2, truth=[]; end
if nargin<3, outputDirectory=''; end
if ischar(source)||isstring(source)
    source=char(source);
    if isfolder(source)
        if isfile(fullfile(source,'simulation_validation.json'))
            report=jsondecode(fileread(fullfile(source,'simulation_validation.json')));
            resultFile=fullfile(report.receiverResultDirectory,'deep_result.mat');
            if isempty(truth), truth=readmatrix(fullfile(source,'truth.csv'),'NumHeaderLines',1); end
        else
            resultFile=fullfile(source,'deep_result.mat');
        end
        if isempty(outputDirectory), outputDirectory=fullfile(source,'navigation_plots'); end
    else
        resultFile=source;
        if isempty(outputDirectory), outputDirectory=fullfile(fileparts(source),'navigation_plots'); end
    end
    result=BDS.B1I.load_result(resultFile);
else
    result=source;
end
t=result.timeSeconds(:);
assert(~isempty(t),'B1I:PlotEmpty','No fusion epochs to plot.');
if ~isempty(outputDirectory)&&~isfolder(outputDirectory), mkdir(outputDirectory); end
figures=gobjects(0); names={};
reference=[];
if ~isempty(truth)
    assert(size(truth,2)==13,'B1I:PlotTruth','Expected 13 truth columns.');
    reference=interp1(truth(:,1),truth(:,8:13),t,'linear',NaN);
end
add('position');
components(result.positionLla,{'Latitude (deg)','Longitude (deg)','Height (m)'},cols(reference,1:3));
add('velocity');
components(result.velocityNed,{'North (m/s)','East (m/s)','Down (m/s)'},cols(reference,4:6));
add('attitude');
components(result.attitudeDeg,{'Yaw (deg)','Pitch (deg)','Roll (deg)'},[]);
add('gyro_bias');
components(result.bias(:,1:3)*180/pi,{'X bias (deg/s)','Y bias (deg/s)','Z bias (deg/s)'},[]);
add('accelerometer_bias');
components(result.bias(:,4:6),{'X bias (m/s^2)','Y bias (m/s^2)','Z bias (m/s^2)'},[]);
add('receiver_clock');
components(result.clockSeconds*1e6,{'Clock bias (us)','Clock drift (us/s)'},[]);
add('trajectory');
if isempty(reference), origin=result.positionLla(1,:); else, origin=reference(1,1:3); end
rotation=dcmecef2ned(origin(1),origin(2)); originEcef=lla2ecef(origin);
ned=zeros(numel(t),3); nedTruth=ned;
for k=1:numel(t)
    ned(k,:)=(lla2ecef(result.positionLla(k,:))-originEcef)*rotation';
    if ~isempty(reference), nedTruth(k,:)=(lla2ecef(reference(k,1:3))-originEcef)*rotation'; end
end
plot(ned(:,2),ned(:,1),'LineWidth',1.2); hold on;
if ~isempty(reference), plot(nedTruth(:,2),nedTruth(:,1),'--','LineWidth',1.2); legend('B1I + IMU','Truth','Location','best'); end
xlabel('East displacement (m)'); ylabel('North displacement (m)'); axis equal; grid on;
title('Common local origin; initial position error retained');
add('tracking_feedback');
subplot(3,1,1); plot(t,result.codeCorrectionChips); ylabel('Code correction (chips)'); grid on;
legend(compose('C%02d',result.prns),'Location','eastoutside');
subplot(3,1,2); plot(t,result.carrierCorrectionHz); ylabel('Carrier correction (Hz)'); grid on;
subplot(3,1,3); stairs(t,sum(result.feedbackApplied,2)); ylabel('Feedback channels'); xlabel('IQ file time (s)'); grid on;
add('tracking_quality');
subplot(3,1,1); plot(t,result.lockMetric); ylabel('Carrier lock metric'); grid on;
legend(compose('C%02d',result.prns),'Location','eastoutside');
subplot(3,1,2); stairs(t,result.measurementCount); ylabel('Accepted satellites'); grid on;
subplot(3,1,3); stairs(t,double(result.measurementUsed)); ylim([-.1,1.1]); ylabel('Filter update (0/1)'); xlabel('IQ file time (s)'); grid on;
if isfield(result,'satnav') && isfield(result,'antennaPositionLla')
    add('antenna_position_comparison');
    components(result.antennaPositionLla,{'Latitude (deg)','Longitude (deg)','Height (m)'}, ...
        result.satnav(:,1:3),{'Fused antenna','Same-loop GNSS'});
    sgtitle('Same-loop diagnostic, not an independent feedback-off control');
    add('antenna_velocity_comparison');
    components(result.antennaVelocityNed,{'North (m/s)','East (m/s)','Down (m/s)'}, ...
        result.satnav(:,4:6),{'Fused antenna','Same-loop GNSS'});
    sgtitle('GNSS diagnostic uses original weighted solver velocity approximation');
end
if isfield(result,'attitudeSigmaDeg')
    add('attitude_standard_deviation');
    components(result.attitudeSigmaDeg,{'Yaw sigma (deg)','Pitch sigma (deg)','Roll sigma (deg)'},[]);
    sgtitle('Filter-reported uncertainty, not measured attitude error');
end
if ~isempty(outputDirectory)
    for k=1:numel(figures)
        exportgraphics(figures(k),fullfile(outputDirectory,[names{k},'.png']),'Resolution',150);
        savefig(figures(k),fullfile(outputDirectory,[names{k},'.fig']));
    end
    fprintf('Saved %d navigation figures (PNG and FIG): %s\n',numel(figures),outputDirectory);
end
    function add(name)
        figures(end+1)=figure('Name',['B1I / INS: ',strrep(name,'_',' ')],'Position',[100 100 1050 750]);
        names{end+1}=name;
    end
    function components(values,labels,expected,legendLabels)
        if nargin<4, legendLabels={'B1I + IMU','Truth'}; end
        for j=1:numel(labels)
            subplot(numel(labels),1,j); plot(t,values(:,j),'LineWidth',1.1); hold on;
            if ~isempty(expected)
                plot(t,expected(:,j),'--','LineWidth',1.1);
                if j==1, legend(legendLabels,'Location','best'); end
            end
            ylabel(labels{j}); grid on;
        end
        xlabel('IQ file time (s)');
    end
end
function values=cols(values,indices)
if ~isempty(values), values=values(:,indices); end
end
