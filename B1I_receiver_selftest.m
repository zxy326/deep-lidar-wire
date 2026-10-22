function B1I_receiver_selftest()
% Synthetic regression suite for MATLAB. No real IQ/IMU or toolbox required.
cfg=B1I_config(); cfg.showPlot=false; cfg.positionAtmosphere=false;
cfg.minElevationDeg=5; cfg.maxBchCorrections=0;
% Carrier/code tracking under 100-Hz acquisition error and 1-ms sign flips.
fixtureFile=[tempname,'.dat'];
cleanFixture=onCleanup(@() delete_fixture(fixtureFile)); %#ok<NASGU>
cfg.iqFile=fixtureFile; cfg.processSeconds=1.1; cfg.startSeconds=0;
samples=0:round(cfg.sampleFreq*0.8)-1; actual=1850; delay=1234;
chip=(samples-delay)*2.046e6*(1+actual/1561.098e6)/cfg.sampleFreq;
primary=BDS.B1I.codeGene(19);
signal=primary(mod(floor(chip),2046)+1).*(-1).^floor(chip/2046);
signal=signal.*exp(2i*pi*(actual*samples/cfg.sampleFreq+0.18));
signal=[signal,zeros(1,round(cfg.sampleFreq*0.3))]; % Deliberate signal loss.
fid=fopen(fixtureFile,'wb','ieee-le'); assert(fid>=0);
closeFixture=onCleanup(@() fclose(fid));
fwrite(fid,[real(signal);imag(signal)],'single'); clear closeFixture;
acq=table(19,delay,1750,'VariableNames',{'PRN','CodeDelaySamples','DopplerHz'});
tracked=BDS.B1I.track(cfg,acq);
assert(numel(tracked.prompt)>1000 && mean(tracked.usable(600:750))>0.9);
assert(abs(median(tracked.dopplerHz(600:750))-actual)<3);
assert(~any(tracked.usable(end-99:end)),'Lost signal was incorrectly marked usable.');
clear signal samples chip;
% Exhaustive systematic BCH messages, strict rejection, optional correction.
for value=0:2047
    message=logical(dec2bin(value,11)-'0'); encoded=bch(message);
    raw=[false(1,15),encoded];
    [out,ok,n]=BDS.B1I.word_decode(raw,true,0);
    assert(ok && n==0 && isequal(out(16:26),message));
    damaged=raw; damaged(16+mod(value,15))=~damaged(16+mod(value,15));
    [~,ok]=BDS.B1I.word_decode(damaged,true,0); assert(~ok);
    [out,ok,n]=BDS.B1I.word_decode(damaged,true,1);
    assert(ok && n==1 && isequal(out,raw));
end

% D1 full symbol stream: NH alignment, inverted polarity, BCH, frame SOW,
% split fields and complete ephemeris assembly. Initial symbol offset = 7 ms.
frames=cell(1,5);
for sf=1:5
    b=header(sf,100002+6*(sf-1));
    if sf==1
        b=put(b,[61,13],1080); b=put(b,[74,9;91,8],12500);
        b=put(b,[226,7;241,17],-13); b=put(b,[258,5;271,17],37);
    elseif sf==2
        b=put(b,[251,12;271,20],round(5153.6*2^19));
        b=put(b,[133,10;151,22],round(0.012*2^33));
        b=put(b,[43,10;61,6],-20); b=put(b,[291,2],floor(12500/2^15));
    elseif sf==3
        b=put(b,[43,10;61,5],mod(12500,2^15));
        b=put(b,[66,17;91,15],round(.94/pi*2^31));
        b=put(b,[212,21;241,11],round(-1.2/pi*2^31));
    end
    frames{sf}=encode_frame(b);
end
bits=[frames{:}]; nh=1-2*('00000100110101001110'-'0')';
prompt=[zeros(7,1);reshape(nh*(1-2*double(bits)),[],1)];
track=fixture(prompt,6,cfg.sampleFreq);
nav=BDS.B1I.decode(track,cfg);
assert(numel(nav.frames)==5 && numel(nav.ephemerides)==1);
e=nav.ephemerides(1);
assert(e.week==1080 && e.toe==100000 && e.toc==100000);
assert(abs(e.af0+13*2^-33)<1e-20 && abs(e.dn+20*pi*2^-43)<1e-20);
assert(abs(e.Omega0+1.2)<2e-9 && nav.symbolOffset==7);
% Streaming refresh: ring-relative epochs must map to the original timeline.
empty=struct('prn',6,'frames',[],'ephemerides',[],'symbolOffset',NaN);
[stream,accepted]=BDS.B1I.refresh_navigation(empty,track,40000,cfg);
assert(accepted && stream.frames(1).epoch==nav.frames(1).epoch+40000);
[again,accepted]=BDS.B1I.refresh_navigation(stream,track,40000,cfg);
assert(accepted && isequaln(stream,again));
[~,accepted]=BDS.B1I.refresh_navigation(stream,track,40001,cfg);
assert(~accepted,'A 1-ms navigation time discontinuity must be rejected.');
% A newly decoded unhealthy flag must revoke even previously valid ephemeris.
healthBits=nav.frames(1).bits; healthBits(43)=true;
healthFrames=frames; healthFrames{1}=encode_frame(healthBits);
healthPrompt=[zeros(7,1);reshape(nh*(1-2*double([healthFrames{:}])),[],1)];
[unhealthy,accepted]=BDS.B1I.refresh_navigation(stream,fixture(healthPrompt,6,cfg.sampleFreq),40000,cfg);
assert(accepted && isempty(unhealthy.ephemerides),'Unhealthy satellite retained valid ephemeris.');

% D2: 10 pages of subframe 1 in 10 successive 3-second frames.
geoFrames=cell(1,50);
for page=1:10
    for sf=1:5
        b=header(sf,200001+3*(page-1));
        if sf==1
            b=put(b,[43,4],page);
            switch page
                case 1
                    b=put(b,[65,13],1080); b=put(b,[78,5;91,12],25000);
                case 3
                    b=put(b,[101,12;121,12],-5);
                    b=put(b,[133,4],15); % af1 upper bits of -1
                case 4
                    b=put(b,[47,6;61,12],2^18-1);
                case 5
                    ecc=round(.001*2^33); b=put(b,[125,10],floor(ecc/2^22));
                case 6
                    ecc=round(.001*2^33); b=put(b,[47,6;61,16],mod(ecc,2^22));
                    b=put(b,[77,6;91,22;121,4],round(6493*2^19));
                case 7
                    b=put(b,[81,2;91,15],25000);
            end
        end
        geoFrames{(page-1)*5+sf}=encode_frame(b);
    end
end
bits=[geoFrames{:}]; prompt=[0;reshape(ones(2,1)*(2*double(bits)-1),[],1)];
navGeo=BDS.B1I.decode(fixture(prompt,1,cfg.sampleFreq),cfg);
assert(numel(navGeo.frames)==10 && numel(navGeo.ephemerides)==1);
g=navGeo.ephemerides(1);
assert(g.toe==200000 && g.af1==-2^-50 && g.isGeo);
[r,dt]=BDS.B1I.satellite(g,g.toe);
assert(abs(norm(r)-g.sqrtA^2*(1-g.ecc))<1e-5 && isfinite(dt));

% A repeated preamble with incorrect word parity must not publish ephemeris.
% Explicit word rejection above covers a full-bit flip; avoid relying on ties
% from flipping only one of the two D2 prompt epochs.
raw=frames{2}; raw(45)=~raw(45);
[~,ok]=BDS.B1I.word_decode(raw(31:60),false,0); assert(~ok);

% Known-position least squares with independently constructed Sagnac ranges.
truth=[6378137,0,0]; clock=830000;
directions=[1,0,0;1,1,0;1,-1,0;1,0,1;1,0,-1;1,1,1];
directions=bsxfun(@rdivide,directions,sqrt(sum(directions.^2,2)));
rs=bsxfun(@plus,truth,directions*22000000); rho=zeros(6,1);
for j=1:6
    theta=7.292115e-5*norm(rs(j,:)-truth)/299792458;
    rot=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1]*rs(j,:)';
    rho(j)=norm(rot'-truth)+clock;
end
fix=BDS.B1I.solve_position(rs,rho,[0,0,0],cfg);
assert(fix.valid && norm(fix.ecef-truth)<0.02 && abs(fix.clockMeters-clock)<0.02);
fix=BDS.B1I.solve_position(rs(1:3,:),rho(1:3),truth,cfg); assert(~fix.valid);
fprintf('B1I receiver synthetic tests passed: BCH, D1, D2, GEO orbit and PVT.\n');
end

function delete_fixture(filename)
if isfile(filename), delete(filename); end
end

function tr=fixture(prompt,prn,fs)
tr=struct('prn',prn,'prompt',complex(prompt),'usable',true(size(prompt)), ...
    'boundarySample',(0:numel(prompt)-1)'*fs/1000);
end

function b=header(sf,sow)
b=false(1,300); b(1:11)=logical('11100010010'-'0');
b=put(b,[16,3],sf); b=put(b,[19,8;31,12],sow);
end

function b=put(b,map,value)
n=sum(map(:,2)); if value<0, value=value+2^n; end
v=logical(dec2bin(value,n)-'0'); offset=0;
for k=1:size(map,1)
    b(map(k,1)+(0:map(k,2)-1))=v(offset+(1:map(k,2))); offset=offset+map(k,2);
end
end

function raw=encode_frame(b)
raw=false(1,300);
for w=1:10
    ix=(w-1)*30+(1:30); data=b(ix);
    if w==1, word=[data(1:15),bch(data(16:26))];
    else, word=reshape([bch(data(1:11));bch(data(12:22))],1,30); end
    raw(ix)=word;
end
end

function encoded=bch(data)
work=[logical(data),false(1,4)];
for k=1:11
    if work(k), work(k+[0,3,4])=~work(k+[0,3,4]); end
end
encoded=[logical(data),work(12:15)];
end
