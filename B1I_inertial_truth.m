function [values,provenance] = B1I_inertial_truth(folder,t)
% Return [YPR deg, gyro bias rad/s, acceleration bias m/s^2].
file=fullfile(folder,'inertial_truth.csv');
if isfile(file)
    data=readmatrix(file,'NumHeaderLines',1);
    assert(size(data,2)==10&&all(isfinite(data(:)))&&all(diff(data(:,1))>0), ...
        'B1I:InertialTruth','Invalid inertial truth file.');
    data(:,2:4)=unwrap(data(:,2:4)*pi/180,[],1)*180/pi;
    values=interp1(data(:,1),data(:,2:10),t,'linear',NaN);
    provenance='Generator-saved inertial_truth.csv';
else
    m=jsondecode(fileread(fullfile(folder,'manifest.json')));
    known='D1 MEO circular orbits; Earth rotation; no atmosphere; FRD body aligned NED; static then accelerate then coast';
    assert(strcmp(m.model,known)&&m.seed==7319&&all(m.config.deep.initialAttitudeDeg==0), ...
        'B1I:MissingInertialTruth','Unknown legacy model: regenerate input with explicit inertial truth.');
    values=repmat([0 0 0,[.003 -.002 .004]*pi/180,.01 -.008 .006],numel(t),1);
    provenance='Legacy fixed-NED generator constants reconstructed from the verified model; not measured truth';
end
assert(all(isfinite(values(:))),'B1I:TruthCoverage','Inertial truth does not cover result times.');
end
