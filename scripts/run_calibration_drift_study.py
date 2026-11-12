"""Study 007: calibration transfer drift. Author: Biswajit Jana."""
import argparse,csv,gzip,hashlib,io,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from calibrated_identification import cp_interval,physical_interval,continuum_test
from calibration_drift import minimax_comparator,signed_drift
from transmission_inference import pt_transmission,wilson

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/study-007';FIG=ROOT/'figures/study-007'
ALLOWANCES=[0.,.002,.004,.006,.008,.01,.015]
BUDGETS=[(100000,30000),(1000000,300000)]

def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))

def write_csv(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def evaluate(y,n,fn,fp,m,k,allowances):
    """Nested acceptance rule; compute binomial intervals once per sample."""
    tail=.05/(6*len(k));pl,pu=cp_interval(y,n,tail);fl,fu=cp_interval(fn,m,tail);gl,gu=cp_interval(fp,m,tail)
    output={};retained=False
    for eps in sorted(set(allowances)):
        if retained:
            output[eps]={'status':'retained_by_nesting','nodes':0};continue
        qlo,qhi=physical_interval(pl,pu,np.maximum(0,fl-eps),np.minimum(1,fu+eps),np.maximum(0,gl-eps),np.minimum(1,gu+eps))
        r=continuum_test(k,qlo,qhi);output[eps]=r
        retained=r['status']!='excluded'
    return output

def case(task):
    index,n,m,d,reps,k,p,q,u=task;seed=2026111700+index;rng=np.random.default_rng(seed)
    methods=[(f'fixed_{eps:g}',eps) for eps in ALLOWANCES]+[('matched_bound',d)]
    counts={(group,method):[0,0,0] for group in ['comparator_null','ellis_power'] for method,eps in methods}
    records=[]
    for group,truth,sign in [('comparator_null',q,1),('ellis_power',p,-1)]:
        fs,gs=signed_drift(.02,u,d,sign);recorded=(1-fs)*truth+gs*(1-truth)
        for rep in range(reps):
            y=rng.binomial(n,recorded);fn=rng.binomial(m,.02,size=24);fp=rng.binomial(m,.02,size=24)
            decisions=evaluate(y,n,fn,fp,m,k,[eps for method,eps in methods])
            row={'group':group,'replicate':rep}
            row.update({f'y{i}':int(v) for i,v in enumerate(y)})
            row.update({f'fn{i}':int(v) for i,v in enumerate(fn)})
            row.update({f'fp{i}':int(v) for i,v in enumerate(fp)})
            for method,eps in methods:
                r=decisions[eps];row[method+'_status']=r['status'];row[method+'_nodes']=r['nodes']
                counts[(group,method)][0]+=r['status']=='excluded'
                counts[(group,method)][1]+=r['status']=='undetermined'
                counts[(group,method)][2]+=r['nodes']
            records.append(row)
    filename=f'case-{index:02d}.csv.gz'
    with (OUT/filename).open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as z:
            with io.TextIOWrapper(z,encoding='utf-8',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
    rates=[]
    for group in ['comparator_null','ellis_power']:
        for method,eps in methods:
            reject,unknown,nodes=counts[(group,method)];lo,hi=wilson(reject,reps)
            rates.append(dict(case=index,N=n,M=m,actual_drift=d,method=method,allowed_drift=eps,
                bound_contains_truth=int(eps>=d),group=group,seed=seed,repetitions=reps,
                rejected=reject,rate=reject/reps,wilson_low=float(lo),wilson_high=float(hi),
                undetermined=unknown,nodes=nodes,file=filename))
    return rates

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight');path=FIG/(name+'.svg')
    fig.savefig(path,bbox_inches='tight',metadata={'Date':None,'Creator':'Biswajit Jana'})
    path.write_bytes(('\n'.join(s.rstrip() for s in path.read_text(encoding='utf-8').splitlines())+'\n').encode());plt.close(fig)

def plots(rates,design):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'study007'})
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
    for col,(n,m) in enumerate(BUDGETS):
        for row,group in enumerate(['comparator_null','ellis_power']):
            ax=axes[row,col]
            for method,label,color in [('fixed_0','Assumes zero drift','#b44f37'),('matched_bound','Allows the scenario drift bound','#126b85')]:
                rr=sorted([r for r in rates if r['N']==n and r['M']==m and r['group']==group and r['method']==method],key=lambda r:r['actual_drift'])
                y=np.array([r['rate'] for r in rr]);lo=np.array([r['wilson_low'] for r in rr]);hi=np.array([r['wilson_high'] for r in rr])
                ax.errorbar([100*r['actual_drift'] for r in rr],y,yerr=[y-lo,hi-y],marker='o',ms=4,capsize=2,label=label,color=color)
            ax.axvline(100*design['matching_drift'],color='gray',ls=':',label='Constructed exchange')
            if row==0:ax.axhline(.05,color='gray',ls='--');ax.set_title(f'N={n:,}; M={m:,} per class/channel')
            if col==0:ax.set_ylabel('Physical-null rejection' if row==0 else 'Ellis rejection power')
            if row==1:ax.set_xlabel('Maximum reference-to-science drift (percentage points)')
            ax.set_ylim(-.02,1.02);ax.grid(alpha=.2);ax.legend(fontsize=8,loc='best')
    fig.suptitle('Structured calibration drift can create or erase apparent discrimination\n500 paired trials per point; pointwise 95% Wilson intervals',fontsize=12)
    fig.tight_layout();save(fig,'01-drift-failure-and-protection')
    fig,axes=plt.subplots(1,2,figsize=(12,5.8));ds=design['drift_grid']
    for ax,group,title in zip(axes,['comparator_null','ellis_power'],['Physical-null rejection','Ellis rejection power']):
        matrix=np.array([[next(r['rate'] for r in rates if r['N']==1000000 and r['group']==group and r['method']==f'fixed_{e:g}' and r['actual_drift']==d) for d in ds] for e in ALLOWANCES])
        im=ax.imshow(matrix,origin='lower',vmin=0,vmax=1,aspect='auto',cmap='Blues')
        ax.set(xticks=range(len(ds)),xticklabels=[f'{100*d:.2f}' for d in ds],yticks=range(len(ALLOWANCES)),yticklabels=[f'{100*e:.1f}' for e in ALLOWANCES],xlabel='Actual drift (percentage points)',ylabel='Assumed maximum drift (percentage points)',title=title)
        ax.tick_params(axis='x',labelrotation=45)
        for j,e in enumerate(ALLOWANCES):
            for i,d in enumerate(ds):
                if e<d:ax.add_patch(Rectangle((i-.5,j-.5),1,1,fill=False,hatch='///',edgecolor='#a6a6a6',linewidth=0))
                ax.text(i,j,f'{100*matrix[j,i]:.0f}',ha='center',va='center',color='white' if matrix[j,i]>.65 else 'black',fontsize=7)
        fig.colorbar(im,ax=ax,label='Rejection fraction')
    fig.suptitle('N=1,000,000; M=300,000 | Hatched cells violate the assumed drift bound\nCell numbers are percentages; all Monte Carlo intervals are in rates.csv',fontsize=12)
    fig.tight_layout();save(fig,'02-assumed-versus-actual-drift')
    k=np.array(design['K']);p=np.array(design['p_ellis']);q=np.array(design['q_comparator']);u=np.array(design['direction'])
    fn,fp=signed_drift(.02,u,design['matching_drift']);recorded=(1-fn)*q+fp*(1-q)
    fig,axes=plt.subplots(1,2,figsize=(12,4.7))
    axes[0].plot(k,.02+.96*p,label='Ellis, unchanged reference response',lw=2)
    axes[0].plot(k,.02+.96*q,'--',label='Comparator, unchanged reference response')
    axes[0].plot(k,recorded,'o',mfc='none',ms=5,label='Comparator after constructed drift')
    axes[0].set(xscale='log',xlabel='K = ka',ylabel='Recorded probability');axes[0].legend(fontsize=8)
    axes[1].plot(k,100*fn,'o-',label='Science false-negative rate');axes[1].plot(k,100*fp,'s-',label='Science false-positive rate')
    axes[1].axhline(2,color='gray',ls='--',label='Reference calibration: 2%')
    axes[1].set(xscale='log',xlabel='K = ka',ylabel='Classification error (%)');axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Population-level construction: identical science and reference-calibration distributions')
    fig.tight_layout();save(fig,'03-exact-drift-exchange')
    # Population geometry, not a Monte Carlo fit or an instrument forecast.
    dstar=design['matching_drift'];grid=np.linspace(0,dstar,101)
    population=[]
    for d in grid:
        population.append(dict(actual_drift=float(d),minimum_allowance_for_overlap=float(dstar-d),
                               matched_allowance=float(d),population_margin_matched=float(dstar-2*d)))
    write_csv(OUT/'population-boundary.csv',population)
    fig,ax=plt.subplots(figsize=(8.2,6))
    x=100*grid;boundary=100*(dstar-grid)
    ax.fill_between(x,boundary,100*dstar,color='#d5e7ef',label='Population overlap with comparator family')
    ax.plot(x,boundary,color='#126b85',lw=2,label='Overlap boundary: allowed + actual drift = d*')
    ax.fill_between(x,0,x,facecolor='none',hatch='///',edgecolor='#b7b7b7',label='Assumed bound smaller than actual drift')
    ax.plot(x,x,'--',color='#b44f37',label='Scenario-matched allowance')
    ax.plot([50*dstar],[50*dstar],'o',color='black')
    ax.annotate(f'Two-model overlap begins at\n{50*dstar:.6f} percentage points',
                (50*dstar,50*dstar),xytext=(.04,.58),arrowprops={'arrowstyle':'->'},fontsize=10)
    ax.set(xlim=(0,100*dstar),ylim=(0,100*dstar),xlabel='Actual Ellis-directed drift (percentage points)',
           ylabel='Allowed drift (percentage points)',title='Infinite-count limit for the prescribed direction\n2% reference errors; saved 24-channel numerical spectra')
    ax.legend(loc='lower left',fontsize=8);ax.grid(alpha=.15)
    fig.tight_layout();save(fig,'04-population-overlap-boundary')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repetitions',type=int,default=500);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    spec=[r for r in read(ROOT/'results/study-002/spectra.csv') if r['design']=='broad_band']
    k=np.array([float(v['K']) for v in spec]);p=np.array([float(v['T_ellis']) for v in spec])
    nominal=json.loads((ROOT/'results/study-002/summary.json').read_text())['design_results']['broad_band']['pt_logC_logb']
    theta,gap,attempts=minimax_comparator(k,p,[nominal,[2.4,.9],[3.,1.1]])
    q=pt_transmission(k,theta);u=(p-q)/gap;dstar=.96*gap
    excluded_gap=gap*(1-1e-5)
    qlo=np.maximum(0,p-excluded_gap-1e-11);qhi=np.minimum(1,p+excluded_gap+1e-11)
    cert=continuum_test(k,qlo,qhi,return_certificate=True)
    assert cert['status']=='excluded',cert
    certificate={'excluded_gap':excluded_gap,'qlo':qlo.tolist(),'qhi':qhi.tolist(),'K':k.tolist(),'p':p.tolist(),**cert}
    (OUT/'population-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8',newline='\n')
    ds=sorted([0.,.001,.002,.003,.004,.005,dstar,.008,.01,.015])
    design={'author':'Biswajit Jana','date':'2026-11-12','repetitions':args.repetitions,'reference_fp_fn':.02,
        'K':k.tolist(),'p_ellis':p.tolist(),'q_comparator':q.tolist(),'theta':theta.tolist(),'direction':u.tolist(),
        'gap':gap,'matching_drift':dstar,'optimization_attempts':attempts,'drift_grid':ds,'allowances':ALLOWANCES,'budgets':BUDGETS,
        'seed_rule':'2026111700 + case index; lower budget first, actual drift increasing',
        'null_pattern':'fn=.02-d*u; fp=.02+d*u','ellis_pattern':'fn=.02+d*u; fp=.02-d*u',
        'selection':'Direction constructed from population spectra, fixed before random counts; adversarial benchmark, not typical instrument drift',
        'alpha':.05,'parameter_domain':{'logC':[-16,25],'b':[.1,30]},
        'source_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['results/study-002/spectra.csv','results/study-002/summary.json']}}
    (OUT/'design.json').write_text(json.dumps(design,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('gap',gap,'matching drift',dstar,'lower gap',excluded_gap,'certificate nodes',cert['nodes'],flush=True)
    tasks=[]
    for n,m in BUDGETS:
        for d in ds:tasks.append((len(tasks),n,m,d,args.repetitions,k,p,q,u))
    rates=[];start=time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(case,t) for t in tasks]
        for i,f in enumerate(as_completed(futures),1):
            r=f.result();rates.extend(r)
            print(f'case {i}/{len(tasks)} N={r[0]["N"]} d={r[0]["actual_drift"]:.6g} {time.time()-start:.1f}s',flush=True)
    rates.sort(key=lambda r:(r['case'],r['group'],r['allowed_drift'],r['method']))
    write_csv(OUT/'rates.csv',rates)
    summary={'author':'Biswajit Jana','experiments':len(tasks)*2*args.repetitions,'count_entries':len(tasks)*2*args.repetitions*72,
        'decision_records':len(rates)*args.repetitions,'undetermined_records':sum(r['undetermined'] for r in rates),
        'max_matched_bound_null_rate':max(r['rate'] for r in rates if r['group']=='comparator_null' and r['method']=='matched_bound'),
        'max_ignored_drift_null_rate':max(r['rate'] for r in rates if r['group']=='comparator_null' and r['method']=='fixed_0'),
        'matching_drift':dstar,'excluded_gap':excluded_gap,
        'data_type':'synthetic counts under a prespecified adversarial drift pattern'}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
    plots(rates,design)
    paths=sorted([*OUT.glob('*.csv'),*OUT.glob('*.gz'),OUT/'design.json',OUT/'summary.json',OUT/'population-certificate.json',*FIG.glob('*.png'),*FIG.glob('*.svg')])
    (OUT/'sha256.json').write_text(json.dumps({str(v.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(v.read_bytes()).hexdigest() for v in paths},indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
