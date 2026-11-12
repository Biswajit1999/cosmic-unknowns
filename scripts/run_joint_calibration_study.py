"""Study 006: joint calibration and full-family exclusion. Author: Biswajit Jana."""
import argparse,csv,gzip,hashlib,io,json,platform,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import scipy,mpmath
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from calibrated_identification import joint_intervals,continuum_test,zero_error_interval
from transmission_inference import pt_transmission,wilson

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/study-006'; FIG=ROOT/'figures/study-006'
SCIENCE=[1000,3000,10000,30000,100000]
CALIBRATION=[0,100,300,1000,3000,10000,30000]

def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))

def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def case(task):
    index,n,m,error,repetitions,k,p,theta=task
    seed=2026111200+index
    rng=np.random.default_rng(seed)
    q=pt_transmission(k,theta)
    records=[];rates=[]
    for group,truth in [('ellis_power',p),('comparator_null',q)]:
        probability=error+(1-2*error)*truth
        rejected=unknown=nodes=0
        for rep in range(repetitions):
            y=rng.binomial(n,probability)
            fn=rng.binomial(m,error,size=24);fp=rng.binomial(m,error,size=24)
            lo,hi=joint_intervals(y,n,fn,fp,m)
            result=continuum_test(k,lo,hi)
            rejected+=result['status']=='excluded'; unknown+=result['status']=='undetermined';nodes+=result['nodes']
            row=dict(group=group,replicate=rep,status=result['status'],nodes=result['nodes'])
            row.update({f'y{i}':int(v) for i,v in enumerate(y)})
            row.update({f'fn{i}':int(v) for i,v in enumerate(fn)})
            row.update({f'fp{i}':int(v) for i,v in enumerate(fp)})
            records.append(row)
        lower,upper=wilson(rejected,repetitions)
        rates.append(dict(N=n,M_per_class_channel=m,true_error=error,group=group,seed=seed,
            rejected=rejected,repetitions=repetitions,rate=rejected/repetitions,
            wilson_low=float(lower),wilson_high=float(upper),undetermined=unknown,total_nodes=nodes,
            science_trials=24*n,calibration_trials=48*m))
    name=f"broad-e{int(error*10000):04d}-N{n:06d}-M{m:06d}.csv.gz"
    with (OUT/name).open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as z:
            with io.TextIOWrapper(z,encoding='utf-8',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
    return rates

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight')
    path=FIG/(name+'.svg');fig.savefig(path,bbox_inches='tight',metadata={'Date':None,'Creator':'Biswajit Jana'})
    path.write_bytes(('\n'.join(s.rstrip() for s in path.read_text(encoding='utf-8').splitlines())+'\n').encode())
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repetitions',type=int,default=500);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();OUT.mkdir(exist_ok=True,parents=True);FIG.mkdir(exist_ok=True,parents=True)
    spec=read(ROOT/'results/study-002/spectra.csv')
    prior=json.loads((ROOT/'results/study-002/summary.json').read_text())
    witnesses=json.loads((ROOT/'results/study-004/summary.json').read_text())
    # Population exclusion brackets; padding represents a numerical sensitivity
    # check, NOT a proved enclosure of the continuum ODE solution.
    population=[];certificates={}
    for design in ['low_band','broad_band']:
        rows=[r for r in spec if r['design']==design]
        k=np.array([float(r['K']) for r in rows]);p=np.array([float(r['T_ellis']) for r in rows])
        for pad in [0.,2e-10]:
            upper=witnesses['designs'][design]['achievable_error_budget']
            lower=upper*.95
            initial=continuum_test(k,*zero_error_interval(p,lower,pad),return_certificate=True)
            assert initial['status']=='excluded'
            best=initial
            for _ in range(10):
                middle=(lower+upper)/2
                r=continuum_test(k,*zero_error_interval(p,middle,pad),return_certificate=True)
                if r['status']=='excluded':lower=middle;best=r
                else:upper=middle
            key=f'{design}-pad{pad:g}'
            certificates[key]={'K':k.tolist(),'p':p.tolist(),'padding':pad,'delta_excluded':lower,
                               'qlo':zero_error_interval(p,lower,pad)[0].tolist(),
                               'qhi':zero_error_interval(p,lower,pad)[1].tolist(),**best}
            population.append(dict(design=design,padding=pad,excluded_through_delta=lower,
                known_feasible_delta=witnesses['designs'][design]['achievable_error_budget'],
                numerical_bracket_upper=upper,nodes=best['nodes']))
            print('population',population[-1],flush=True)
    write_csv(OUT/'population-brackets.csv',population)
    (OUT/'exclusion-certificates.json').write_text(json.dumps(certificates,indent=2)+'\n',encoding='utf-8',newline='\n')
    rows=[r for r in spec if r['design']=='broad_band']
    k=np.array([float(r['K']) for r in rows]);p=np.array([float(r['T_ellis']) for r in rows])
    theta=np.array(prior['design_results']['broad_band']['pt_logC_logb'])
    tasks=[]
    for error in [0.,.01]:
        for n in SCIENCE:
            for m in CALIBRATION:tasks.append((len(tasks),n,m,error,args.repetitions,k,p,theta))
    design={'author':'Biswajit Jana','date':'2026-11-12','data_type':'synthetic science and known-class calibration counts',
            'alpha':.05,'tail_probability':.05/144,'science_budgets':SCIENCE,'calibration_budgets':CALIBRATION,
            'repetitions_per_group':args.repetitions,'true_fp_fn':[0.,.01],'K':k.tolist(),'p_ellis':p.tolist(),
            'null_theta':theta.tolist(),'seed_rule':'2026111200 + case index; zero error first, N then M increasing',
            'parameter_domain':{'logC':[-16,25],'b':[.1,30]},'max_nodes':12000,
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mpmath.__version__},
            'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['results/study-002/spectra.csv','results/study-002/summary.json','results/study-004/summary.json']}}
    (OUT/'design.json').write_text(json.dumps(design,indent=2)+'\n',encoding='utf-8',newline='\n')
    rates=[];started=time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(case,t) for t in tasks]
        for i,future in enumerate(as_completed(futures),1):
            result=future.result();rates.extend(result)
            print(f"case {i}/{len(tasks)}: N={result[0]['N']} M={result[0]['M_per_class_channel']} e={result[0]['true_error']} power={result[0]['rate']:.3f}; {time.time()-started:.1f}s",flush=True)
    rates.sort(key=lambda r:(r['true_error'],r['N'],r['M_per_class_channel'],r['group']))
    write_csv(OUT/'rates.csv',rates)
    summary={'author':'Biswajit Jana','experiments':len(tasks)*2*args.repetitions,
        'science_channel_counts':len(tasks)*2*args.repetitions*24,
        'calibration_channel_counts':len(tasks)*2*args.repetitions*48,
        'undetermined':sum(r['undetermined'] for r in rates),
        'null_rejections':sum(r['rejected'] for r in rates if r['group']=='comparator_null'),
        'null_experiments':len(tasks)*args.repetitions,
        'interpretation':'Joint-confidence-region exclusion of a bounded comparator family; not physical wormhole existence.'}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'study006'})
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),sharey=True)
    for ax,error in zip(axes,[0.,.01]):
        for m in [0,100,300,1000,3000,30000]:
            rr=sorted([r for r in rates if r['true_error']==error and r['M_per_class_channel']==m and r['group']=='ellis_power'],key=lambda r:r['N'])
            y=np.array([r['rate'] for r in rr]);lo=np.array([r['wilson_low'] for r in rr]);hi=np.array([r['wilson_high'] for r in rr])
            ax.errorbar([r['N'] for r in rr],y,yerr=[y-lo,hi-y],marker='o',ms=3,capsize=2,label=f'M={m:,}')
        ax.set(xscale='log',ylim=(-.02,1.02),xlabel='Science quanta per channel, N',title=f'True false-positive = false-negative = {100*error:g}%')
        ax.grid(alpha=.2);ax.legend(title='Calibration / class / channel',fontsize=8)
    axes[0].set_ylabel('Full-family rejection fraction')
    fig.suptitle(f'Independent calibration can restore discrimination | {args.repetitions} trials/point; pointwise 95% Wilson bars')
    fig.tight_layout();save(fig,'01-calibration-power')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for ax,error in zip(axes,[0.,.01]):
        matrix=np.array([[next(r['rate'] for r in rates if r['true_error']==error and r['N']==n and r['M_per_class_channel']==m and r['group']=='ellis_power') for n in SCIENCE] for m in CALIBRATION])
        im=ax.imshow(matrix,origin='lower',aspect='auto',vmin=0,vmax=1,cmap='Blues')
        ax.set(xticks=range(len(SCIENCE)),xticklabels=[f'{n:,}' for n in SCIENCE],yticks=range(len(CALIBRATION)),yticklabels=[f'{m:,}' for m in CALIBRATION],xlabel='N per science channel',ylabel='M per calibration class/channel',title=f'True error per class: {100*error:g}%')
        for j in range(len(CALIBRATION)):
            for i in range(len(SCIENCE)):ax.text(i,j,f'{100*matrix[j,i]:.1f}%',ha='center',va='center',color='white' if matrix[j,i]>.6 else 'black',fontsize=8)
        fig.colorbar(im,ax=ax,label='Power')
    fig.suptitle('Power depends on science and calibration budgets; cells are discrete simulated designs')
    fig.tight_layout();save(fig,'02-budget-map')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=True)
    for ax,error in zip(axes,[0.,.01]):
        rr=[r for r in rates if r['true_error']==error and r['group']=='comparator_null']
        y=np.array([r['rate'] for r in rr]);lo=np.array([r['wilson_low'] for r in rr]);hi=np.array([r['wilson_high'] for r in rr])
        ax.errorbar(range(len(rr)),y,yerr=[y-lo,hi-y],fmt='o',ms=3,capsize=2)
        ax.axhline(.05,color='gray',ls='--',label='Theoretical size upper bound')
        ax.set(title=f'True error per class: {100*error:g}%',xlabel='Grid case (N then M increasing)',ylim=(-.002,.06));ax.legend(fontsize=8);ax.grid(alpha=.2)
    axes[0].set_ylabel('Comparator-null rejection fraction')
    fig.suptitle('Null checks at one physical comparator; uniform claim follows from coverage, not this sample')
    fig.tight_layout();save(fig,'03-null-check')
    files=sorted([*OUT.glob('*.csv'),*OUT.glob('*.gz'),OUT/'summary.json',OUT/'design.json',OUT/'exclusion-certificates.json',*FIG.glob('*.png'),*FIG.glob('*.svg')])
    hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (OUT/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
