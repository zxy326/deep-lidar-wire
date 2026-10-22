function reports=B1I_regression_suite(dataRoot)
% One entry for component checks, optionally followed by all four IQ fixtures.
% No argument: component tests only. With dataRoot: expensive full replays.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
B1I_selftest();
B1I_code_variance_selftest();
B1I_source_quality_selftest();
B1I_consistency_selftest();
B1I_source_contract_selftest();
B1I_bootstrap_prior_selftest();
B1I_partial_update_selftest();
B1I_lifecycle_selftest();
B1I_dynamic_week_selftest();
B1I_weak_integration_selftest();
B1I_deep_selftest();
B1I_receiver_selftest();
reports={};
if nargin==0
    fprintf('Component tests passed. Full IQ replay/real-data acceptance NOT run.\n');
    return;
end
cases={'turn_7320','turn_7321','vertical_7320','vertical_7321'};
for k=1:numel(cases)
    folder=fullfile(dataRoot,cases{k});
    assert(isfile(fullfile(folder,'manifest.json')),'B1I:FixtureMissing','Missing fixture %s.',folder);
end
for k=1:numel(cases)
    fprintf('Full regression %d/%d: %s\n',k,numel(cases),cases{k});
    reports{k}=B1I_consistency_replay(fullfile(dataRoot,cases{k})); %#ok<AGROW>
end
fprintf('All four structural replay regressions completed; inspect reported accuracy changes separately.\n');
end
