function B1I_weight_validation(seed)
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
for scenario={'turn','vertical'}
 folder=fullfile(root,'~temp','b1i_repeat_20261004',sprintf('%s_%d',scenario{1},seed));
 report=B1I_weight_trial(folder);
 fid=fopen(fullfile(folder,'source_weight_report.json'),'w'); assert(fid>=0);
 fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true)); fclose(fid); close all;
end
end
