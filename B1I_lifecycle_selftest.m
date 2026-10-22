function B1I_lifecycle_selftest()
root=fileparts(mfilename('fullpath')); addpath(genpath(root));
h=BDS.B1I.PromptHistory(.003);
for k=1:5, h.push(k,k*10,mod(k,2)==1); end
[snapshot,offset]=h.snapshot(19);
assert(offset==2 && isequal(snapshot.prompt,[3;4;5]) && isequal(snapshot.boundarySample,[30;40;50]));
cfg=B1I_deep_config(); cfg.sampleFreq=4092000; cfg.startSeconds=0; cfg.processSeconds=3.1;
cfg.deep.navigationWindowSeconds=1; cfg.deep.weakHoldSeconds=.5;
cfg.iqFile=[tempname,'.dat']; fileCleanup=onCleanup(@()remove_file(cfg.iqFile));
fs=cfg.sampleFreq; n=round(fs*cfg.processSeconds); ix=0:n-1; code=BDS.B1I.codeGene(19);
frequency=1850; delay=1234; chips=(ix-delay)*2.046e6*(1+frequency/1561.098e6)/fs;
x=code(mod(floor(chips),2046)+1).*(-1).^floor(chips/2046).*exp(2i*pi*(frequency*ix/fs+.18));
% Brief blackout then a real signal return; later blackout exceeds hold limit.
x(ix/fs>=.85 & ix/fs<1.05)=0; x(ix/fs>=1.6 & ix/fs<2.4)=0;
fid=fopen(cfg.iqFile,'wb','ieee-le'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fwrite(fid,[real(x);imag(x)],'single'); clear closer x ix chips;
acq=table(19,delay,1850,'VariableNames',{'PRN','CodeDelaySamples','DopplerHz'});
ch=BDS.B1I.DeepChannel(cfg,acq); channelCleanup=onCleanup(@()delete(ch));
ch.advance(.7); assert(ch.locked); ch.beginFusion(); ch.feedback(0,0);
everHeld=false;
for t=.71:.01:1.5
    ch.advance(t);
    if ch.holdActive
        everHeld=true;
        [~,~,valid]=ch.observe(t,struct('sow',100000,'epoch',1));
        assert(~valid,'A held channel must never masquerade as a GNSS measurement.');
    end
    if ~ch.stopped, ch.feedback(0,0); end
end
assert(everHeld && ~ch.stopped && ch.locked,'Brief blackout did not recover.');
for t=1.51:.01:2.4
    ch.advance(t); if ~ch.stopped, ch.feedback(0,0); end
end
assert(ch.stopped && ~ch.locked,'Continued aiding must not defeat the finite loss timeout.');
% Exercise the same acquisition/replacement function used by the receiver.
blank=BDS.B1I.reacquire_channel(ch,2.3,cfg.processSeconds,cfg);
assert(isempty(blank),'A blank IQ window must not acquire a satellite.');
scan=cfg; scan.dopplerMaxHz=2000; scan.dopplerStepHz=250;
replacement=BDS.B1I.reacquire_channel(ch,2.5,cfg.processSeconds,scan);
assert(~isempty(replacement)); replacementCleanup=onCleanup(@()delete(replacement));
assert(replacement.ptr/cfg.sampleFreq<=2.5 && ~replacement.navigationConfirmed);
% New-PRN discovery uses the same past-only acquisition, independent of old.
[newCandidates,newScan]=BDS.B1I.past_acquisition(scan,[19,20],2.5,cfg.processSeconds);
assert(any(newCandidates.PRN==19),'Dynamic search did not discover present PRN.');
assert(newScan.startSeconds+.002*cfg.noncoherentPeriods<=2.5);
replacement.advance(3.05);
assert(replacement.locked && ~replacement.stopped,'Returned signal did not reacquire.');
[~,~,valid]=replacement.observe(3.05,struct('sow',100000,'epoch',1));
assert(~valid,'Reacquisition must wait for fresh navigation before publishing measurements.');
% Stale acceleration assistance must expire even if the channel stays locked.
replacement.beginFusion(); replacement.carrierRate=10;
replacement.lastFeedbackTime=2.4; replacement.advance(3.06);
assert(replacement.carrierRate==0,'Expired aiding acceleration was retained.');
replacement.cfg.deep=rmfield(replacement.cfg.deep,'weakHoldSeconds');
replacement.advance(3.07); % Old saved configurations use the safe age default.
% Unvalidated reacquisition must not publish observables, even when locked.
ch.stopped=false; ch.locked=true; ch.navigationConfirmed=false;
[~,~,valid]=ch.observe(ch.ptr/fs,struct('sow',100000,'epoch',1)); assert(~valid);
clear replacementCleanup channelCleanup fileCleanup;
fprintf('B1I history, brief-loss recovery, bounded hold, IQ reacquisition and navigation quarantine tests passed.\n');
end
function remove_file(file)
if isfile(file), delete(file); end
end
