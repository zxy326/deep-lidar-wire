function fig=B1I_plot_weak_tracking(source)
% Inspect recorded estimator decisions, not configured signal truth.
if ischar(source)||isstring(source)
    if isfolder(source), source=fullfile(source,'deep_result.mat'); end
    r=BDS.B1I.load_result(source);
else
    r=source;
end
assert(isfield(r,'coherentMilliseconds') && isfield(r,'cn0DbHz'), ...
    'B1I:OldResult','Rerun with adaptive-integration diagnostics enabled.');
fig=figure('Name','B1I weak-signal integration');
tiledlayout(3,1);
nexttile; plot(r.timeSeconds,r.cn0DbHz); grid on;
ylabel('Estimated C/N0 (dB-Hz)'); yline(37,'--'); yline(39,':');
legend(compose('C%02d',r.prns),'Location','eastoutside');
nexttile; stairs(r.timeSeconds,r.coherentMilliseconds); grid on;
ylabel('Coherent time (ms)'); yticks([1 2 4]);
nexttile; stairs(r.timeSeconds,sum(r.nhSynchronized,2)); grid on;
ylabel('Verified symbol/NH channels'); xlabel('IQ file time (s)');
end
