function B1I_dynamic_week_selftest()
% Fast non-IQ tests. Lifecycle test separately exercises actual acquisition.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
e=struct('week',1081,'sow',604770,'toe',604000,'toc',604000);
assert(BDS.B1I.ephemeris_age(e,604801,1081)==801);
assert(BDS.B1I.ephemeris_age(e,1,1082)==801);
assert(BDS.B1I.ephemeris_age(e,1,1083)>604800,'Old-week alias accepted.');
e.sow=0; e.week=1082; e.toe=604792; e.toc=604792;
assert(BDS.B1I.ephemeris_age(e,604801,1081)==9,'Previous-week toe incorrectly rejected.');
e.week=0; e.toe=0; e.toc=0;
assert(BDS.B1I.ephemeris_age(e,604801,8191)==1,'13-bit week carry failed.');
% Wrapping transmitted SOW must not jump the pseudorange/time difference.
r=[604799.99,604800,604800.01]; raw=mod(r-.08,604800);
tx=r+mod(raw-r+302400,604800)-302400;
assert(max(abs((r-tx)-.08))<2e-10);
% Explicit storage schema guards non-channel Nx6/Nx9 state arrays.
s=struct('timeSeconds',(1:4)','prns',[19,20],'channelAddedTimeSeconds',[0,0]);
for n={'feedbackApplied','feedbackEligible','feedbackNonzero','nhSynchronized'}
    s.(n{1})=true(4,2);
end
for n={'codeCorrectionChips','carrierCorrectionHz','lockMetric','coherentMilliseconds','cn0DbHz'}
    s.(n{1})=ones(4,2);
end
s.diagnostics=struct('biasAfter',ones(4,6),'updateReason',zeros(4,1,'uint8'));
for n={'rangeGateResidualMeters','rateGateResidualMps','postRangeResidualMeters', ...
        'postRateResidualMps','channelBoundaryLagSeconds','rawPseudorangeMeters', ...
        'filterPseudorangeMeters','rangeVarianceMeters2','rateVarianceMps2'}
    s.diagnostics.(n{1})=ones(4,2);
end
for n={'accepted','used','rangeAccepted','rateAccepted','weakHold'}
    s.diagnostics.(n{1})=true(4,2);
end
s.diagnostics.reason=zeros(4,2,'uint8'); before=s;
s=BDS.B1I.append_channel_storage(s,21,2.5);
assert(isequal(s.prns,[19,20,21]) && s.channelAddedTimeSeconds(3)==2.5);
assert(isequal(s.diagnostics.biasAfter,before.diagnostics.biasAfter));
assert(isequal(s.lockMetric(:,1:2),before.lockMetric) && all(isnan(s.lockMetric(:,3))));
assert(all(~s.feedbackApplied(:,3)) && all(s.diagnostics.reason(:,3)==1));
rejected=false;
try, BDS.B1I.append_channel_storage(s,21,3); catch err, rejected=strcmp(err.identifier,'B1I:DuplicatePRN'); end
assert(rejected);
fprintf('B1I continuous week/ephemeris age and dynamic column tests passed.\n');
end
