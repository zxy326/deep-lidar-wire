function stats=B1I_height_report(dataDirectory)
% Plot/measure height from a saved simulation run; never inject height truth.
summary=jsondecode(fileread(fullfile(dataDirectory,'simulation_validation.json')));
r=BDS.B1I.load_result(fullfile(summary.receiverResultDirectory,'deep_result.mat'));
truth=readmatrix(fullfile(dataDirectory,'truth.csv'),'NumHeaderLines',1);
height=interp1(truth(:,1),truth(:,10),r.timeSeconds,'linear',NaN);
errorHeight=r.positionLla(:,3)-height;
assert(~isempty(errorHeight)&&all(isfinite(errorHeight)),'B1I:HeightTruth','Invalid height truth coverage.');
stats=struct('heightRmseMeters',sqrt(mean(errorHeight.^2)), ...
    'maxHeightErrorMeters',max(abs(errorHeight)),'initialHeightErrorMeters',errorHeight(1), ...
    'finalHeightErrorMeters',errorHeight(end),'epochs',numel(errorHeight), ...
    'firstTimeSeconds',r.timeSeconds(1),'lastTimeSeconds',r.timeSeconds(end), ...
    'updatedEpochs',nnz(r.measurementUsed),'unlockedChannelEpochs',nnz(r.diagnostics.reason==4), ...
    'minimumAcceptedSatellites',min(r.measurementCount),'prns',r.prns, ...
    'receiverResultDirectory',summary.receiverResultDirectory);
f=figure('Name','B1I height truth validation','Position',[100 100 1000 650]);
subplot(2,1,1); plot(r.timeSeconds,r.positionLla(:,3),'LineWidth',1.2); hold on;
plot(r.timeSeconds,height,'--','LineWidth',1.2); ylabel('Height (m)'); legend('B1I + IMU','Truth'); grid on;
subplot(2,1,2); plot(r.timeSeconds,errorHeight,'LineWidth',1.2); yline(0,'k--');
ylabel('Height error (m)'); xlabel('IQ file seconds'); grid on;
exportgraphics(f,fullfile(dataDirectory,'height_validation.png'));
savefig(f,fullfile(dataDirectory,'height_validation.fig'));
fid=fopen(fullfile(dataDirectory,'height_validation.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',jsonencode(stats,'PrettyPrint',true));
fprintf('Height RMSE %.4f m, maximum %.4f m, final %.4f m; %d/%d updates\n', ...
    stats.heightRmseMeters,stats.maxHeightErrorMeters,stats.finalHeightErrorMeters,stats.updatedEpochs,stats.epochs);
end
