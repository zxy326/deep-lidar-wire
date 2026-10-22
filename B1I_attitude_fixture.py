"""Paired rotating-body IMU fixture, reusing zero-lever-arm antenna IQ."""
from pathlib import Path
import argparse
import json
import numpy as np
from B1I_simulate import body_motion, make_imu

def check_model():
    # Independent numerical matrix derivative validates Euler/body-rate signs.
    for t in np.linspace(24.01,80,101):
        _,r,rate=body_motion(t)
        h=1e-5
        derivative=(body_motion(t+h)[1]-body_motion(t-h)[1])/(2*h)
        skew=r.T@derivative
        numerical=np.array([skew[2,1],skew[0,2],skew[1,0]])
        assert np.linalg.norm(rate-numerical)<1e-8
        assert np.linalg.norm(r.T@r-np.eye(3))<1e-12
    assert np.all(body_motion(23)[0]==0)

def generate(base,output):
    check_model()
    base=Path(base).resolve(); output=Path(output).resolve()
    m=json.loads((base/'manifest.json').read_text(encoding='utf-8'))
    assert np.linalg.norm(m['config']['deep']['armBodyMeters'])==0
    assert m['scenario']=='straight', 'This paired trial expects the straight baseline.'
    assert not output.exists(), 'Refusing to overwrite a fixture.'
    duration=m['config']['processSeconds']
    imu,truth=make_imu(duration,np.random.default_rng(m['seed']),m['scenario'],True)
    original=np.loadtxt(base/'truth.csv',delimiter=',',skiprows=1)
    assert np.max(np.abs(truth-original))<1e-7, 'Antenna trajectory changed.'
    output.mkdir(parents=True)
    np.savetxt(output/'imu.csv',imu,delimiter=',',fmt='%.12f')
    (output/'truth.csv').write_bytes((base/'truth.csv').read_bytes())
    inertial=np.column_stack((truth[:,0],np.array([body_motion(t)[0] for t in truth[:,0]]),
        np.tile(np.deg2rad([.003,-.002,.004]),(len(truth),1)),
        np.tile([.01,-.008,.006],(len(truth),1))))
    np.savetxt(output/'inertial_truth.csv',inertial,delimiter=',',fmt='%.12f',
        header='iq_seconds,yaw_deg,pitch_deg,roll_deg,bgx_rad_s,bgy_rad_s,bgz_rad_s,bax_m_s2,bay_m_s2,baz_m_s2',comments='')
    m['config']['imuFile']=str(output/'imu.csv')
    m['config']['outputDir']=str(output/'matlab_runs')
    m['attitudeMotion']='Smooth ZYX yaw/pitch/roll after 24 s; zero lever arm; same translational truth and IQ'
    m['pairedBaseDirectory']=str(base)
    m['model']=m['model'].replace('FRD body aligned NED','FRD body with analytic ZYX rotation')
    (output/'manifest.json').write_text(json.dumps(m,indent=2),encoding='utf-8')
    print('Attitude derivative tests passed; paired IMU generated:',output)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();generate(a.base,a.output)
