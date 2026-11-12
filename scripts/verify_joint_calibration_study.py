"""Independent artifact checks for study 006. Author: Biswajit Jana."""
import csv,gzip,hashlib,json
from pathlib import Path
import numpy as np
from mpmath import iv
from calibrated_identification import joint_intervals,continuum_test
from transmission_inference import pt_transmission,wilson

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/study-006'

def main():
    design=json.loads((OUT/'design.json').read_text())
    for name,digest in json.loads((OUT/'sha256.json').read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    for name,digest in design['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    # Recheck every population leaf at twice the production precision.
    iv.dps=70
    certs=json.loads((OUT/'exclusion-certificates.json').read_text())
    leaf_total=0
    for certificate in certs.values():
        k=certificate['K'];lo=certificate['qlo'];hi=certificate['qhi']
        D=iv.mpf(certificate['delta_excluded']);pad=iv.mpf(certificate['padding'])
        # Check that the saved, padded physical bounds enclose the exact
        # population response inversion, not merely the subsequent search.
        for target,lbound,ubound in zip(certificate['p'],lo,hi):
            P=iv.mpf(target)
            exact_lower=(P-pad-D)/(1-D);exact_upper=(P+pad)/(1-D)
            if bool(exact_lower.a>0):assert bool((iv.mpf(lbound)-exact_lower).b<=0)
            else:assert lbound==0
            if bool(exact_upper.b<1):assert bool((iv.mpf(ubound)-exact_upper).a>=0)
            else:assert ubound==1
        leaves=sorted(certificate['leaves'])
        assert leaves[0][0]==.1 and leaves[-1][1]==30.
        for previous,current in zip(leaves,leaves[1:]):assert previous[1]==current[0]
        for bl,br,i,j in leaves:
            def value(b,index,q):
                z=iv.mpf(float(b))*iv.mpf(float(k[index]))
                v=iv.mpf(float(q))
                return 2*iv.log((iv.exp(z)-iv.exp(-z))/2)-iv.log(v/(1-v))
            l=iv.mpf(-16) if i==0 else value(bl,i-1,hi[i-1])
            u=iv.mpf(25) if j==0 else value(br,j-1,lo[j-1])
            assert bool((l-u).a>0)
            leaf_total+=1
    witness=json.loads((ROOT/'results/study-004/summary.json').read_text())
    upper_certificates={}
    for band in ['low_band','broad_band']:
        cert=certs[f'{band}-pad0'];w=witness['designs'][band]
        # Small explicit slack converts a floating-point construction into a
        # checked feasible allowance for the fixed saved target probabilities.
        delta=w['achievable_error_budget']+1e-10
        theta=w['theta_logC_logb'];C=iv.exp(iv.mpf(theta[0]));b=iv.exp(iv.mpf(theta[1]));D=iv.mpf(delta)
        for kvalue,pvalue in zip(cert['K'],cert['p']):
            z=b*iv.mpf(kvalue);s=(iv.exp(z)-iv.exp(-z))/2;qiv=s*s/(s*s+C);piv=iv.mpf(pvalue)
            assert bool(((1-D)*qiv-piv).b<=0)
            assert bool((qiv+D*(1-qiv)-piv).a>=0)
        upper_certificates[band]={'verified_feasible_delta':delta,'precision_digits':70,
                                  'target':'exact saved binary floating-point probabilities'}
    iv.dps=35
    with (OUT/'rates.csv').open() as f:rates=list(csv.DictReader(f))
    k=np.array(design['K']);p=np.array(design['p_ellis']);q=pt_transmission(k,design['null_theta'])
    total=audited=0;idx=0
    for error in design['true_fp_fn']:
        for n in design['science_budgets']:
            for m in design['calibration_budgets']:
                seed=2026111200+idx;idx+=1;rng=np.random.default_rng(seed)
                name=f'broad-e{int(error*10000):04d}-N{n:06d}-M{m:06d}.csv.gz'
                with gzip.open(OUT/name,'rt') as f:records=list(csv.DictReader(f))
                for group,truth in [('ellis_power',p),('comparator_null',q)]:
                    subset=[r for r in records if r['group']==group]
                    assert len(subset)==design['repetitions_per_group']
                    rejected=unknown=nodes=0
                    for rep,r in enumerate(subset):
                        y=np.array([int(r[f'y{i}']) for i in range(24)])
                        fn=np.array([int(r[f'fn{i}']) for i in range(24)])
                        fp=np.array([int(r[f'fp{i}']) for i in range(24)])
                        assert np.array_equal(y,rng.binomial(n,error+(1-2*error)*truth))
                        assert np.array_equal(fn,rng.binomial(m,error,size=24))
                        assert np.array_equal(fp,rng.binomial(m,error,size=24))
                        assert int(r['replicate'])==rep
                        rejected+=r['status']=='excluded';unknown+=r['status']=='undetermined';nodes+=int(r['nodes']);total+=1
                        # Deterministic sample of fresh solver replays across every case/group.
                        if rep%50==0:
                            result=continuum_test(k,*joint_intervals(y,n,fn,fp,m))
                            assert result['status']==r['status'] and result['nodes']==int(r['nodes'])
                            audited+=1
                    rate=next(v for v in rates if float(v['true_error'])==error and int(v['N'])==n and int(v['M_per_class_channel'])==m and v['group']==group)
                    assert rejected==int(rate['rejected']) and unknown==int(rate['undetermined']) and nodes==int(rate['total_nodes'])
                    np.testing.assert_allclose(wilson(rejected,len(subset)),[float(rate['wilson_low']),float(rate['wilson_high'])],atol=1e-14)
    summary=json.loads((OUT/'summary.json').read_text())
    assert total==summary['experiments']
    result={'experiments_replayed':total,'count_entries_replayed':total*72,'solver_decisions_replayed':audited,
            'population_certificate_leaves_verified_at_70_digits':leaf_total,'artifact_hashes':'all verified',
            'upper_witness_certificates':upper_certificates}
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))

if __name__=='__main__':main()
