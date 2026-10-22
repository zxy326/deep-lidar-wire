function report=B1I_weak_integration_selftest(seed)
% NH/boundary tests plus sampled IQ comparison. No trajectory truth in receiver.
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
if nargin==0, seed=7361; end
rng(seed);
nh=1-2*('00000100110101001110'-'0');
for prn=[19,1]
    for anchor=1:20
        h=BDS.B1I.CoherentIntegrator(prn); h.synchronize(anchor);
        h.requested=min(4,h.period); completed=0;
        for e=anchor:anchor+79
            slot=mod(e-anchor,h.period); bit=(-1)^floor((e-anchor)/h.period);
            signNH=1; if h.period==20, signNH=nh(slot+1); end
            [ready,z,T,df]=h.push(e,bit*signNH*[2 1 .5],.001,3,true);
            if ready
                assert(abs(abs(z(1))-2*h.length)<1e-12);
                assert(abs(T-h.length*.001)<1e-12 && abs(df-3*T)<1e-12);
                assert(mod(slot+1,h.length)==0); completed=completed+1;
            end
        end
        assert(completed==80/h.length);
    end
end
h=BDS.B1I.CoherentIntegrator(19); h.requested=4;
[ready,~,T]=h.push(5,[1 1 1],.001,0,true);
assert(ready && T==.001 && h.length==1,'Unverified NH must never enable long integration.');
h.synchronize(1); h.requested=4; h.push(1,[1 1 1],.001,0,true);
h.synchronize(21); assert(h.count==1,'Same-phase navigation refresh must preserve partial window.');
h.synchronize(22); assert(h.count==0 && h.length==1,'Changed anchor must flush partial integration.');
h=BDS.B1I.CoherentIntegrator(19); h.synchronize(1); h.requested=4;
h.push(1,[1 1 1],.001,0,true); h.clearWindow();
for e=2:4
    ready=h.push(e,[1 1 1],.001,0,true);
    assert(~ready,'An interrupted coherent window must not publish a partial measurement.');
end
for e=5:8
    [ready,z]=h.push(e,nh(e)*[1 1 1],.001,0,true);
end
assert(ready && abs(z(1)-4)<1e-12,'Feedback-boundary restart failed.');
% Moment estimator and hysteresis use noisy prompts, not supplied C/N0 truth.
h=BDS.B1I.CoherentIntegrator(19); h.synchronize(1);
for e=1:2000
    p=sqrt(10^(34/10)*.001)+complex(randn,randn)/sqrt(2);
    h.push(e,[p p p],.001,0,true);
end
assert(h.length==4 && abs(h.cn0-34)<4,'Weak C/N0 decision failed.');
for e=2001:2400
    p=sqrt(10^(48/10)*.001)+complex(randn,randn)/sqrt(2);
    h.push(e,[p p p],.001,0,true);
end
assert(h.length==1 && h.cn0>39,'Strong signal did not return to 1 ms.');

cfg=B1I_deep_config(); cfg.sampleFreq=4092000; cfg.startSeconds=0; cfg.processSeconds=2.4;
cfg.deep.navigationWindowSeconds=0; cfg.iqFile=[tempname,'.dat'];
cleanup=onCleanup(@()remove_file(cfg.iqFile));
fs=cfg.sampleFreq; frequency=500; delay=1234; code=BDS.B1I.codeGene(19);
fid=fopen(cfg.iqFile,'wb','ieee-le'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
for k=0:2399
    ix=k*4092+(0:4091); chips=(ix-delay)*2.046e6*(1+frequency/1561.098e6)/fs;
    epoch=floor(chips/2046); bits=(-1).^floor(epoch/20);
    signNH=nh(mod(epoch,20)+1);
    x=code(mod(floor(chips),2046)+1).*bits.*signNH.*exp(2i*pi*frequency*ix/fs);
    cn0=48; if k>=800 && k<2000, cn0=34; end
    noise=sqrt(fs/10^(cn0/10)/2);
    x=x+noise*(randn(size(x))+1i*randn(size(x)));
    fwrite(fid,[real(x);imag(x)],'single');
end
clear closer;
acq=table(19,delay,frequency,'VariableNames',{'PRN','CodeDelaySamples','DopplerHz'});
report=struct('seed',seed,'cn0StrongDbHz',48,'cn0WeakDbHz',34, ...
    'note','Single-channel synchronized IQ with ideal frequency aid; not end-to-end positioning acceptance.');
for adaptive=[false,true]
    cfg.deep.adaptiveIntegration=adaptive;
    ch=BDS.B1I.DeepChannel(cfg,acq); channelCleanup=onCleanup(@()delete(ch));
    ch.advance(.7); assert(ch.locked); ch.coherent.synchronize(1); ch.beginFusion();
    long=0; valid=0; total=0;
    for t=.71:.01:2.4
        ch.advance(t);
        if ch.stopped, break; end
        % Ideal bounded inertial frequency aid, identical for both receivers.
        ch.feedback(0,(frequency-ch.carrInt)*299792458/1561.098e6);
        long=long+(ch.coherent.length>1); valid=valid+ch.locked; total=total+1;
    end
    result=struct('stopped',ch.stopped,'lockedEpochs',valid,'processedEpochs',total, ...
        'longIntegrationEpochs',long,'finalCN0',ch.coherent.cn0,'finalMs',ch.coherent.length);
    if adaptive, report.adaptive=result; else, report.fixed=result; end
    clear channelCleanup;
end
assert(~report.adaptive.stopped && report.adaptive.longIntegrationEpochs>50, ...
    'Adaptive channel did not sustain the weak IQ interval.');
assert(report.adaptive.finalMs==1,'Recovery did not return to short integration.');
assert(report.adaptive.lockedEpochs>report.fixed.lockedEpochs, ...
    'This seeded weak-signal comparison did not demonstrate a lock benefit.');
out=fullfile(root,'~temp','b1i_weak_integration'); if ~isfolder(out), mkdir(out); end
fid=fopen(fullfile(out,sprintf('seed_%d.json',seed)),'w'); assert(fid>=0);
closer=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fprintf('%s\n',jsonencode(report));
fprintf('B1I NH boundaries, adaptive integration and sampled weak IQ tests passed.\n');
end
function remove_file(path)
if isfile(path), delete(path); end
end
