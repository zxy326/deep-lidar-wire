function results=B1I_repeat_validation(parent,seeds)
% Fixed-parameter paired comparisons, retaining each successful/failed case.
if nargin<1, parent=fullfile(fileparts(mfilename('fullpath')),'~temp','b1i_repeat_20261004'); end
if nargin<2, seeds=[7320,7321]; end
summaryName=['repeat_summary',sprintf('_%d',seeds),'.json'];
results=struct([]); n=0;
for scenario={'turn','vertical'}
 for seed=seeds
  n=n+1; folder=fullfile(parent,sprintf('%s_%d',scenario{1},seed));
  cache=fullfile(folder,'paired_report.json');
  row=struct('scenario',scenario{1},'seed',seed,'error','','report',struct());
  try
   assert(isfile(fullfile(folder,'manifest.json')),'B1I:Fixture','Missing generated fixture.');
   if isfile(cache)
    row.report=jsondecode(fileread(cache)); fprintf('Using saved paired report: %s\n',cache);
   else
    row.report=B1I_feedback_trial(folder);
    writejson(cache,row.report);
   end
  catch err
   row.error=getReport(err,'extended','hyperlinks','off'); warning('%s',row.error);
  end
  if n==1, results=row; else, results(end+1)=row; end %#ok<AGROW>
  writejson(fullfile(parent,summaryName),results);
  close all;
 end
end
end
function writejson(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
