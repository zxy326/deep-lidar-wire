function reports=B1I_repeat_report(parent,outputDirectory)
% Read all four saved trials; no replay, receiver interpolation or exclusions.
if nargin<1, parent=fullfile(fileparts(mfilename('fullpath')),'~temp','b1i_repeat_20261004'); end
if nargin<2, outputDirectory=parent; end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
f=figure('Name','Repeated feedback comparisons','Position',[100 100 1200 800]);
reports=cell(4,1); n=0;
for scenario={'turn','vertical'}
 for seed=[7320,7321]
  n=n+1; folder=fullfile(parent,sprintf('%s_%d',scenario{1},seed));
  r=jsondecode(fileread(fullfile(folder,'paired_report.json')));
  reports{n}=struct('scenario',scenario{1},'seed',seed,'report',r);
  s=load(fullfile(r.outputDirectory,'feedback_ab_inputs.mat'),'B','C');
  truth=readmatrix(fullfile(folder,'truth.csv'),'NumHeaderLines',1);
  subplot(2,2,n); hold on;
  runs={s.B,s.C};
  for j=1:2
   run=runs{j}; t=run.timeSeconds(:); ref=interp1(truth(:,1),truth(:,5:7),t);
   v=zeros(numel(t),3);
   for k=1:numel(t)
    v(k,:)=run.velocityNed(k,:)*dcmecef2ned(run.positionLla(k,1),run.positionLla(k,2));
   end
   plot(t,vecnorm(v-ref,2,2));
  end
  title(sprintf('%s, seed %d',scenario{1},seed)); xlabel('IQ seconds'); ylabel('3D velocity error (m/s)');
  grid on; legend('No feedback','Acceleration-aided feedback');
  fprintf('%s %d: P %.6f -> %.6f; H %.6f -> %.6f; V %.6f -> %.6f; completed %d/%d\n', ...
   scenario{1},seed,r.B.positionRmseMeters,r.C.positionRmseMeters,r.B.heightRmseMeters,r.C.heightRmseMeters, ...
   r.B.velocityRmseMps,r.C.velocityRmseMps,r.B.completed,r.C.completed);
 end
end
exportgraphics(f,fullfile(outputDirectory,'B1I_repeat_velocity.png'));
savefig(f,fullfile(outputDirectory,'B1I_repeat_velocity.fig'));
fid=fopen(fullfile(outputDirectory,'B1I_repeat_results.json'),'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(reports,'PrettyPrint',true));
end
