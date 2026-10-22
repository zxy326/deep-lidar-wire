"""Generate synchronized synthetic B1I D1 fc32 IQ, IMU and independent truth.

Requires NumPy. Synthetic circular broadcast orbits, not a real almanac.
No receiver positions/observations are injected into the tested receiver.
"""
from pathlib import Path
import argparse
import json
import hashlib
import re
import shutil
import numpy as np

C=299792458.; RF=1561098000.; RATE=2046000.; WE=7.292115e-5
MU=3.986004418e14; A=6378137.; E2=6.6943799901413165e-3
TOE=90000.; START=89999.; NAVSTART=89994.

def ecef(lat,lon,h):
    lat,lon=np.deg2rad([lat,lon]); n=A/np.sqrt(1-E2*np.sin(lat)**2)
    return np.array([(n+h)*np.cos(lat)*np.cos(lon),(n+h)*np.cos(lat)*np.sin(lon),
                     (n*(1-E2)+h)*np.sin(lat)])

def lla(r):
    p=np.hypot(r[0],r[1]); lat=np.arctan2(r[2],p*(1-E2))
    for _ in range(12):
        n=A/np.sqrt(1-E2*np.sin(lat)**2); lat=np.arctan2(r[2]+E2*n*np.sin(lat),p)
    n=A/np.sqrt(1-E2*np.sin(lat)**2)
    return np.array([lat,np.arctan2(r[1],r[0]),p*np.cos(lat)+r[2]*np.sin(lat)-n*(1-E2*np.sin(lat)**2)])

def ned(lat,lon):
    s,c=np.sin(lat),np.cos(lat); sl,cl=np.sin(lon),np.cos(lon)
    return np.array([[-s*cl,-s*sl,c],[-sl,cl,0],[-c*cl,-c*sl,-s]])

R0=ecef(30,115,100); C0=ned(*np.deg2rad([30,115]))

def trajectory(t,scenario='straight'):
    # Static until 24 s; constant horizontal acceleration for 8 s; then coast.
    t=np.asarray(t); u=np.clip(t-24,0,8); after=np.maximum(t-32,0)
    accel=np.array([.4,.25,0.])
    displacement=(.5*u*u+8*after)[...,None]*accel
    velocity=u[...,None]*accel
    acceleration=((t>=24)&(t<32))[...,None]*accel
    if scenario!='straight':
        u=np.maximum(t-24,0); w=np.pi/12; active=(t>=24)
        if scenario=='turn':
            # Turning translation with body fixed to NED (not a yaw manoeuvre).
            displacement=np.stack([8*(1-np.cos(w*u)),8*(w*u-np.sin(w*u)),0*u],axis=-1)
            velocity=np.stack([8*w*np.sin(w*u),8*w*(1-np.cos(w*u)),0*u],axis=-1)
            acceleration=np.stack([8*w*w*np.cos(w*u),8*w*w*np.sin(w*u),0*u],axis=-1)*active[...,None]
        elif scenario=='vertical':
            w=2*np.pi/12
            displacement=np.stack([.1*u*u,0*u,-2*(1-np.cos(w*u))],axis=-1)
            velocity=np.stack([.2*u,0*u,-2*w*np.sin(w*u)],axis=-1)
            acceleration=np.stack([.2+0*u,0*u,-2*w*w*np.cos(w*u)],axis=-1)*active[...,None]
        else: raise ValueError('Unknown trajectory scenario: '+scenario)
    return R0+displacement@C0,velocity@C0,acceleration@C0

def orbit(e,t):
    dt=np.asarray(t)-TOE; a=e['sqrtA']**2
    u=e['M0']+np.sqrt(MU/a**3)*dt
    node=e['Omega0']+(e['OmegaDot']-WE)*dt-WE*TOE
    x=a*np.cos(u); y=a*np.sin(u); ci=np.cos(e['i0']); si=np.sin(e['i0'])
    return np.stack([x*np.cos(node)-y*ci*np.sin(node),x*np.sin(node)+y*ci*np.cos(node),y*si],axis=-1)

def range_model(e,t,scenario='straight'):
    receiver,_,_=trajectory(t,scenario); tau=np.full(np.asarray(t).shape,.08)
    for _ in range(6):
        satellite=orbit(e,START+np.asarray(t)-tau); th=WE*tau
        rotated=np.stack([np.cos(th)*satellite[...,0]+np.sin(th)*satellite[...,1],
                          -np.sin(th)*satellite[...,0]+np.cos(th)*satellite[...,1],satellite[...,2]],axis=-1)
        distance=np.linalg.norm(rotated-receiver,axis=-1); tau=distance/C
    return distance,rotated

def constellation():
    candidates=[]
    for node in np.linspace(-np.pi,np.pi,36,endpoint=False):
        for phase in np.linspace(-np.pi,np.pi,36,endpoint=False):
            e=dict(sqrtA=round(5282.6*2**19)/2**19,i0=round(np.deg2rad(55)/np.pi*2**31)*np.pi/2**31,
                   Omega0=round(node/np.pi*2**31)*np.pi/2**31,M0=round(phase/np.pi*2**31)*np.pi/2**31,
                   OmegaDot=round(-8e-9/np.pi*2**43)*np.pi/2**43)
            distance,r=range_model(e,0.); los=(r-R0)/distance; elev=np.arcsin(-los@C0[2])
            if elev>np.deg2rad(15): candidates.append((e,los,elev))
    # Greedy D-optimal selection, beginning with a high-elevation satellite.
    selected=[max(range(len(candidates)),key=lambda j:candidates[j][2])]
    while len(selected)<6:
        def score(j):
            H=np.array([np.r_[-candidates[k][1],1] for k in selected+[j]])
            return np.linalg.slogdet(H.T@H+np.eye(4)*1e-7)[1]
        selected.append(max((j for j in range(len(candidates)) if j not in selected),key=score))
    out=[]
    for prn,j in zip(range(19,25),selected):
        e=candidates[j][0].copy(); e['prn']=prn; e['elevationDeg']=float(np.rad2deg(candidates[j][2])); out.append(e)
    H=np.array([np.r_[-candidates[j][1],1] for j in selected])
    return out,float(np.sqrt(np.trace(np.linalg.inv(H.T@H))))

def put(bits,mapping,value):
    width=sum(n for _,n in mapping); value=int(value)
    assert -(1<<(width-1))<=value<(1<<width)
    value%=1<<width; encoded=np.array(list(f'{value:0{width}b}'),dtype=int); k=0
    for start,n in mapping: bits[start-1:start-1+n]=encoded[k:k+n]; k+=n

def bch(data):
    work=np.r_[data,np.zeros(4,int)]
    for k in range(11):
        if work[k]:work[k+np.array([0,3,4])]^=1
    return np.r_[data,work[-4:]]

def frame(e,sow):
    sf=(sow//6)%5+1; b=np.zeros(300,int); b[:11]=np.array(list('11100010010'),int)
    put(b,[(16,3)],sf);put(b,[(19,8),(31,12)],sow % 604800)
    if sf==1:
        put(b,[(61,13)],(1081+sow//604800) % 8192);put(b,[(74,9),(91,8)],int(TOE/8))
    elif sf==2:
        put(b,[(93,20),(121,12)],round(e['M0']/np.pi*2**31))
        put(b,[(251,12),(271,20)],round(e['sqrtA']*2**19))
        put(b,[(291,2)],int(TOE/8)//2**15)
    elif sf==3:
        put(b,[(43,10),(61,5)],int(TOE/8)%2**15)
        put(b,[(66,17),(91,15)],round(e['i0']/np.pi*2**31))
        put(b,[(132,11),(151,13)],round(e['OmegaDot']/np.pi*2**43))
        put(b,[(212,21),(241,11)],round(e['Omega0']/np.pi*2**31))
    raw=[]
    for w in range(10):
        d=b[w*30:(w+1)*30]
        raw.extend(np.r_[d[:15],bch(d[15:26])] if w==0 else np.column_stack([bch(d[:11]),bch(d[11:22])]).ravel())
    return np.array(raw,int)

def primary(prn,repo):
    source=(repo/'+BDS/+B1I/codeGene.m').read_text(encoding='ascii')
    matrix=source.split('taps = [ ...')[1].split('];')[0]
    taps=[[int(x) for x in row.split()] for row in matrix.strip().split(';')][prn-1]
    r1=r2=int('01010101010',2)
    m1=sum(1<<x for x in [0,1,2,3,4,10]);m2=sum(1<<x for x in [0,2,3,6,7,8,9,10])
    phase=sum(1<<(11-x) for x in taps if x) # stage s maps to bit (11-s)
    result=[]
    for _ in range(2046):
        result.append(2*((r1^((r2&phase).bit_count()&1))&1)-1)
        r1=(r1>>1)|(((r1&m1).bit_count()&1)<<10);r2=(r2>>1)|(((r2&m2).bit_count()&1)<<10)
    return np.array(result)

def body_motion(t):
    """Smooth ZYX attitude after 24 s; R maps body-column vectors to NED."""
    u=max(float(t)-24.,0.); w=2*np.pi/20
    amp=np.deg2rad([20.,8.,6.]); freq=w*np.array([1.,1.3,.7])
    yaw,pitch,roll=amp*(1-np.cos(freq*u))
    yd,pd,rd=amp*freq*np.sin(freq*u)
    cy,sy=np.cos(yaw),np.sin(yaw); cp,sp=np.cos(pitch),np.sin(pitch)
    cr,sr=np.cos(roll),np.sin(roll)
    R=np.array([[cy*cp,cy*sp*sr-sy*cr,cy*sp*cr+sy*sr],
                [sy*cp,sy*sp*sr+cy*cr,sy*sp*cr-cy*sr],[-sp,cp*sr,cp*cr]])
    rate=np.array([rd-yd*sp,pd*cr+yd*sr*cp,-pd*sr+yd*cr*cp])
    return np.rad2deg([yaw,pitch,roll]),R,rate

def make_imu(duration,rng,scenario='straight',attitude_motion=False):
    times=np.arange(round(duration*100)+1)/100
    positions,velocities,accelerations=trajectory(times,scenario); imu=[]; truth=[]
    for t,r,v,a in zip(times,positions,velocities,accelerations):
        lat,lon,h=lla(r); cn=ned(lat,lon); vn=v@cn.T; an=a@cn.T
        n=A/np.sqrt(1-E2*np.sin(lat)**2); m=A*(1-E2)/(1-E2*np.sin(lat)**2)**1.5
        earth=WE*np.array([np.cos(lat),0,-np.sin(lat)])
        transport=np.array([vn[1]/(n+h),-vn[0]/(m+h),-vn[1]*np.tan(lat)/(n+h)])
        gravity=9.7803253359*(1+.00193185265241*np.sin(lat)**2)/np.sqrt(1-E2*np.sin(lat)**2)*(1-2*h/A)
        gyro=earth+transport
        force=an+np.cross(2*earth,vn)-np.array([0,0,gravity])
        if attitude_motion:
            _,rotation,body_rate=body_motion(t)
            gyro=rotation.T@gyro+body_rate
            force=rotation.T@force
        gyro+=np.deg2rad([.003,-.002,.004])+rng.normal(0,np.deg2rad(.002),3)
        force+=np.array([.01,-.008,.006])+rng.normal(0,.003,3)
        imu.append(np.r_[1000+t,gyro,force]); truth.append(np.r_[t,r,v,np.rad2deg([lat,lon]),h,vn])
    return np.array(imu),np.array(truth)

def filtered_code_tables(code,nh,fs,oversample=8):
    """Continuous-delay, zero-phase bandlimited rectangular-code model.

    Fourier coefficients include the analytic rectangular chip pulse, not
    point-sampled floor(chips). Three 20-ms symbols provide edge context;
    return the centre symbol for four possible neighbouring polarities.
    This is an idealized front end, not a measured hardware filter.
    """
    base=np.tile(code,20)*np.repeat(nh,2046)
    n=len(base)*3; length=n*oversample
    k=np.arange(length//2+1); frequency=k*RATE/n
    stop=min(.90*fs/2,.95*RATE); start=.80*stop
    response=np.zeros_like(frequency)
    response[frequency<=start]=1
    taper=(frequency>start)&(frequency<stop)
    response[taper]=.5*(1+np.cos(np.pi*(frequency[taper]-start)/(stop-start)))
    pulse=np.sinc(k/n)*np.exp(-1j*np.pi*k/n)*response*oversample
    tables=[]
    for prev,nxt in ((1,1),(1,-1),(-1,1),(-1,-1)):
        symbols=np.concatenate((prev*base,base,nxt*base))
        spectrum=np.fft.fft(symbols)
        wave=np.fft.irfft(spectrum[k%n]*pulse,n=length)
        centre=wave[len(base)*oversample:2*len(base)*oversample+1]
        tables.append(centre)
    tables=np.array(tables)
    # Unit mean signal power, so the stated C/N0 convention is retained.
    tables/=np.sqrt(np.mean(tables[:,:-1]**2))
    return tables

def filtered_spread(transmit,bits,tables,oversample=8):
    symbol=np.floor(transmit*50).astype(np.int64)
    current=2*bits[symbol]-1
    prev=(2*bits[symbol-1]-1)*current
    nxt=(2*bits[symbol+1]-1)*current
    row=2*(prev<0).astype(int)+(nxt<0).astype(int)
    phase=np.mod(transmit*RATE,40920)*oversample
    index=np.floor(phase).astype(np.int64); fraction=phase-index
    return current*((1-fraction)*tables[row,index]+fraction*tables[row,index+1])

def generate(output,repo,duration=36.,fs=4092000,seed=7319,signal_model='bandlimited',scenario='straight',late_prn_seconds=0.):
    output=Path(output).resolve(); output.mkdir(parents=True,exist_ok=True)
    iqfile=output/'b1i_sim_fc32.dat'
    if iqfile.exists(): raise FileExistsError(f'Refusing to overwrite {iqfile}')
    expected=round(duration*fs)*8
    if shutil.disk_usage(output).free<expected+500_000_000:raise RuntimeError('Insufficient free space')
    rng=np.random.default_rng(seed); eph,gdop=constellation()
    imu,truth=make_imu(duration,rng,scenario)
    np.savetxt(output/'imu.csv',imu,delimiter=',',fmt='%.12f')
    np.savetxt(output/'truth.csv',truth,delimiter=',',fmt='%.12f',header='iq_seconds,x,y,z,vx,vy,vz,latitude,longitude,height,vn,ve,vd',comments='')
    # Explicit independent attitude/bias truth, never passed to the receiver.
    inertial_truth=np.column_stack((truth[:,0],np.zeros((len(truth),3)),
        np.tile(np.deg2rad([.003,-.002,.004]),(len(truth),1)),
        np.tile([.01,-.008,.006],(len(truth),1))))
    np.savetxt(output/'inertial_truth.csv',inertial_truth,delimiter=',',fmt='%.12f',
        header='iq_seconds,yaw_deg,pitch_deg,roll_deg,bgx_rad_s,bgy_rad_s,bgz_rad_s,bax_m_s2,bay_m_s2,baz_m_s2',comments='')
    bits=[np.concatenate([frame(e,s) for s in range(int(NAVSTART),int(START+duration)+12,6)]) for e in eph]
    codes=[primary(e['prn'],repo) for e in eph];nh=1-2*np.array(list('00000100110101001110'),int)
    assert signal_model in ('bandlimited','legacy')
    tables=[filtered_code_tables(code,nh,fs) for code in codes] if signal_model=='bandlimited' else None
    cn0=52.;sigma=np.sqrt(fs/(2*10**(cn0/10))); knots=np.arange(round(duration*100)+1)/100
    ranges=[range_model(e,knots,scenario)[0] for e in eph]
    chunk=round(fs*.01); digest=hashlib.sha256()
    print(f'Generating {duration}s, {fs} Hz, {expected/1e9:.3f} GB, 6 satellites, GDOP={gdop:.2f}',flush=True)
    with iqfile.open('xb') as f:
        for start in range(0,round(duration*fs),chunk):
            t=(start+np.arange(min(chunk,round(duration*fs)-start)))/fs
            signal=np.zeros(t.size,complex)
            for j,e in enumerate(eph):
                distance=np.interp(t,knots,ranges[j]); transmit=t+(START-NAVSTART)-distance/C
                chips=transmit*RATE; ms=np.floor(transmit*1000).astype(np.int64)
                symbols=2*bits[j][np.floor(transmit*50).astype(np.int64)]-1
                if signal_model=='bandlimited':
                    spread=filtered_spread(transmit,bits[j],tables[j])
                else:
                    spread=codes[j][np.floor(chips).astype(np.int64)%2046]*nh[ms%20]*symbols
                if j==len(eph)-1 and late_prn_seconds>0:
                    spread=spread*(t>=late_prn_seconds)
                signal+=spread*np.exp(2j*np.pi*(-distance/(C/RF)+.13*j))
            signal+=sigma*(rng.normal(size=t.size)+1j*rng.normal(size=t.size))
            raw=signal.astype('<c8').tobytes();f.write(raw);digest.update(raw)
            if start%round(fs*5)==0:print(f'IQ {start/fs:.0f}/{duration:.0f}s',flush=True)
    cfg=dict(iqFile=str(iqfile),imuFile=str(output/'imu.csv'),sampleFreq=fs,processSeconds=duration,
             codeReplica='bandlimited' if signal_model=='bandlimited' else 'rectangular',
             codeReplicaStopHz=min(.90*fs/2,.95*RATE),
             startSeconds=0,prns=[e['prn'] for e in eph],maxChannels=6,showPlot=False,
             acquisitionScanOffsets=[0],dopplerMaxHz=6000,positionAtmosphere=False,
             outputDir=str(output/'matlab_runs'),
             deep=dict(bootstrapSeconds=24,dt=.01,initialAttitudeDeg=[0,0,0],armBodyMeters=[0,0,0]),
             imu=dict(format='csv',gyroUnits='rad/s',accelUnits='m/s^2',specificForceConfirmed=True),
             sync=dict(confirmed=True,imuReferenceSeconds=1000,imuToIqOffsetSeconds=0,imuClockScale=1,imuAxesToBody=np.eye(3).tolist()))
    manifest=dict(config=cfg,seed=seed,scenario=scenario,cn0DbHz=cn0,receiverBdtStart=START,bdtWeek=1081,ephemerides=eph,
                  latePrn=eph[-1]['prn'],latePrnSeconds=late_prn_seconds,
                  signalModel=signal_model,signalModelVersion=2 if signal_model=='bandlimited' else 1,
                  frontEnd=dict(description='Analytic rectangular chip spectrum, zero-phase raised-cosine lowpass, neighbour-symbol context',
                                stopHz=min(.90*fs/2,.95*RATE),lookupSamplesPerChip=8) if tables is not None else None,
                  iqSha256=digest.hexdigest(),iqBytes=expected,gdopAtStart=gdop,
                  truthColumns=['iq_seconds','x','y','z','vx','vy','vz','latitude','longitude','height','vn','ve','vd'],
                  model='D1 MEO circular orbits; Earth rotation; no atmosphere; FRD body aligned NED; trajectory='+scenario,
                  validation='Generated input only. Not proof MATLAB receiver passed.')
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Generated synchronized IQ / IMU / truth:',output,flush=True)
    return manifest

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);parser.add_argument('--seconds',type=float,default=36)
    parser.add_argument('--signal-model',choices=['bandlimited','legacy'],default='bandlimited')
    parser.add_argument('--seed',type=int,default=7319)
    parser.add_argument('--scenario',choices=['straight','turn','vertical'],default='straight')
    parser.add_argument('--late-prn-seconds',type=float,default=0.)
    parser.add_argument('--week-rollover',action='store_true')
    args=parser.parse_args()
    if args.week_rollover: TOE,START,NAVSTART=604768.,604769.,604764.
    generate(args.output,Path(__file__).resolve().parent,args.seconds,seed=args.seed,signal_model=args.signal_model,scenario=args.scenario,late_prn_seconds=args.late_prn_seconds)
