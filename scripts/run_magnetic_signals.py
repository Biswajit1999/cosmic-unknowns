"""Finite-band and delayed-pulse study. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_hermitenorm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from magnetic_signals import (mouth_probability, amplitudes, matching_solve,
    local_transmission, halfwidth, eventual_fraction, collected_fraction,
    collection_pulses, gaussian_average, gaussian_fourier, gaussian_lorentzian, pulse_amplitude,
    compact_envelope, compact_autocorrelation)
from magnetic_energy import MagneticEnergy, scale_diagnostics

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/study-014'; FIG=ROOT/'figures/study-014'


def write_csv(name, rows):
    with (OUT/name).open('w', encoding='utf-8', newline='') as stream:
        writer=csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def save(fig, name):
    fig.savefig(FIG/(name+'.png'), dpi=180, bbox_inches='tight', metadata={'Author':'Biswajit Jana'})
    p=FIG/(name+'.svg'); fig.savefig(p, bbox_inches='tight', metadata={'Date':None,'Creator':'Biswajit Jana'})
    p.write_text('\n'.join(l.rstrip() for l in p.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8', newline='\n')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'Biswajit-Jana-study014'})
    spectra=[]; matching=[]
    for x in [.005,.02,.05]:
        for det in np.linspace(-.5,.5,1001):
            spectra.append(dict(x=x, detuning=float(det), transmission=float(abs(local_transmission(x,det))**2)))
        for z in [.25,.5,1.,1.001,1.01,1.5,2.,2.123]:
            r,t=amplitudes(x,z); rr,tt=matching_solve(x,z)
            matching.append(dict(x=x,z=z,flux_error=float(abs(abs(r)**2+abs(t)**2-1)),
                matching_error=float(max(abs(r-rr),abs(t-tt)))))
    write_csv('spectra.csv',spectra); write_csv('matching_audit.csv',matching)
    bands=[]; band_audit=[]
    for x in [.005,.02,.05]:
        for b in np.geomspace(halfwidth(x)*.001,10,101):
            v,err=gaussian_average(x,float(b))
            bands.append(dict(x=x,bandwidth=float(b),bandwidth_over_hwhm=float(b/halfwidth(x)),
                average=v,quad_error_estimate=err,lorentzian=float(gaussian_lorentzian(x,b)),
                duration_sigma_over_ell=float(1/(2*b))))
        for b in [.005,.02,.1,.5,2.]:
            for det in [0.,.2,.5]:
                v,err=gaussian_average(x,b,det); independent=gaussian_fourier(x,b,det)
                band_audit.append(dict(x=x,bandwidth=b,detuning=det,quadrature=v,fourier=independent,
                    absolute_difference=abs(v-independent),quad_error_estimate=err))
    write_csv('bandwidth.csv',bands); write_csv('bandwidth_audit.csv',band_audit)
    collection=[]; delays=[]
    for x in np.geomspace(1e-6,.1,121):
        p=float(mouth_probability(x)); n=float(collection_pulses(x))
        delays.append(dict(x=float(x),mouth_probability=p,hwhm_z=float(halfwidth(x)),first_pulse=p*p,
            eventual=float(eventual_fraction(x)),n90_continuous=n,n90_complete=int(np.ceil(n)),
            t90_over_ell=float((2*np.ceil(n)-1)*np.pi),
            sigma_t_over_ell_at_bandwidth_hwhm=float(1/(2*halfwidth(x)))))
    for x in [.005,.02,.05]:
        for n in np.unique(np.ceil(np.geomspace(1,1e6,401))).astype(int):
            total=float(collected_fraction(x,n))
            collection.append(dict(x=x,pulses=int(n),arrival_over_pi_ell=int(2*n-1),
                incident_energy_fraction=total,fraction_of_eventual=total/float(eventual_fraction(x))))
    write_csv('collection.csv',collection); write_csv('delay_scaling.csv',delays)
    u=np.linspace(0,10*np.pi,10001); amplitude=pulse_amplitude(u,.05,10.,pulses=6)
    pulses=[dict(u_over_pi=float(v/np.pi),energy_flux=float(abs(a)**2),flux_over_first_energy=float(abs(a)**2/mouth_probability(.05)**2)) for v,a in zip(u,amplitude)]
    write_csv('pulse_profile.csv',pulses)
    # Direct coherent time integrals, including pulse overlap; compare a finite
    # frequency-domain geometric series with the exact Fourier overlap sum.
    time_checks=[]
    for b in [.15,.5,2.]:
        x=.05; count=40; p=float(mouth_probability(x)); R=1-p
        domain=(-10/b,(2*count-1)*np.pi+10/b)
        centers=[(2*m+1)*np.pi for m in range(count)]
        integral,error=quad(lambda v: float(abs(pulse_amplitude(np.array(v),x,b,count))**2),
            *domain,points=centers,epsabs=1e-12,epsrel=1e-10,limit=1000)
        m=np.arange(count); differences=m[:,None]-m[None,:]
        overlap=p*p*np.sum(R**(m[:,None]+m[None,:])*np.exp(-2*np.pi**2*b*b*differences**2))
        time_checks.append(dict(bandwidth=b,pulses=count,time_integral=integral,overlap_sum=float(overlap),
            absolute_difference=float(abs(integral-overlap)),quad_error_estimate=error))
    write_csv('time_audit.csv',time_checks)
    compact=[]
    for duration in [.4,1.4,8.]:
        for chirp in [0.,7.]:
            x=.05;p=float(mouth_probability(x));R=1-p;count=12
            centers=(2*np.arange(count)+1)*np.pi
            def flux(u):
                amp=-p*np.sum(R**np.arange(count)*compact_envelope(u-centers,duration,chirp))
                return float(abs(amp)**2)
            points=sorted(set(np.concatenate([centers-duration/2,centers,centers+duration/2])))
            direct=sum(quad(flux,l,r,epsabs=1e-13,epsrel=1e-11)[0] for l,r in zip(points[:-1],points[1:]))
            baseline=float(collected_fraction(x,count));prediction=baseline
            for lag in range(1,count):
                corr=compact_autocorrelation(2*np.pi*lag,duration,chirp)
                prediction+=2*p*p*R**lag*(1-R**(2*(count-lag)))/(1-R*R)*corr.real
            compact.append(dict(duration_over_ell=duration,chirp=chirp,pulses=count,
                envelope_norm=compact_autocorrelation(0,duration,chirp).real,
                first_roundtrip_autocorrelation=compact_autocorrelation(2*np.pi,duration,chirp).real,
                direct_energy=direct,autocorrelation_prediction=prediction,separated_prediction=baseline,
                absolute_difference=abs(direct-prediction)))
    write_csv('compact_pulse_audit.csv',compact)
    # Real parameter algebra from Study 013, but this is not a certified solution.
    eq=MagneticEnergy().minimum(.01); scenario=[]; physical=[]
    for q in np.geomspace(1e6,1e28,111):
        scales=scale_diagnostics(q,.3,.01,eq['y'],eq['energy'])
        binding_re=scales['binding_times_ell']/scales['ell_over_re']
        x=.1*binding_re; n=float(collection_pulses(x)); t=(2*np.ceil(n)-1)*np.pi
        physical.append(dict(q=float(q),x=x,binding_energy_times_re=binding_re,
            ell_over_re=scales['ell_over_re'],d_over_ell=np.pi/eq['y'],
            orbital_omega_ell=scales['orbital_omega_ell'],omega_times_t90=scales['orbital_omega_ell']*t,
            t90_over_ell=t,bandwidth_for_sigma_t_equal_d=eq['y']/(2*np.pi),
            resonance_hwhm=float(halfwidth(x))))
    for q in [1e9,1e15,1e21,1e27]:
        scales=scale_diagnostics(q,.3,.01,eq['y'],eq['energy']); L=scales['ell_over_re']
        binding_re=scales['binding_times_ell']/L; x=.1*binding_re; b=1000.
        nodes,weights=roots_hermitenorm(128); weights/=np.sqrt(2*np.pi)
        xx=x+b*nodes/L
        assert min(xx)>0 and max(xx)<.2
        averaged_first=float(weights@(mouth_probability(xx)**2))
        averaged_eventual=float(weights@eventual_fraction(xx))
        n90=int(np.ceil(collection_pulses(x)))
        averaged_collected=float(weights@collected_fraction(xx,n90))
        scenario.append(dict(q=q,g=.3,s=.01,x=x,bandwidth_z=b,ell_over_re=L,
            binding_energy_times_re=binding_re,quantum_energy_over_binding=.1,
            sigma_t_over_d=1/(2*b)/(np.pi/eq['y']),bandwidth_over_carrier=b/(x*L),
            eventual_frozen=float(eventual_fraction(x)),eventual_variable_x=averaged_eventual,
            relative_eventual_change=averaged_eventual/float(eventual_fraction(x))-1,
            first_frozen=float(mouth_probability(x)**2),first_variable_x=averaged_first,
            relative_first_change=averaged_first/float(mouth_probability(x)**2)-1,
            collected_at_n90_variable_x=averaged_collected,
            relative_collected_change=averaged_collected/float(collected_fraction(x,n90))-1,
            sigma_t_at_hwhm_over_d=1/(2*halfwidth(x))/(np.pi/eq['y']),
            n90=float(collection_pulses(x)),omega_times_t90=scales['orbital_omega_ell']*(2*np.ceil(collection_pulses(x))-1)*np.pi))
    write_csv('parameter_scan.csv',physical); write_csv('scenario_audit.csv',scenario)
    # Figure 1: resolve peaks locally; coarse comb is illustrative only.
    fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for x,c in zip([.005,.02,.05],['#166b8a','#b16535','#55906c']):
        rr=[r for r in bands if r['x']==x]
        ax[0].loglog([r['bandwidth_over_hwhm'] for r in rr],[r['average'] for r in rr],color=c,label=f'omega r_e = {x:g}')
        ax[0].axhline(eventual_fraction(x),color=c,ls=':',lw=.8)
    ax[0].set(xlabel='Energy bandwidth / resonance half-width',ylabel='Total transmitted / incident radial energy',title='A perfect peak does not imply a perfect packet')
    ax[0].legend(fontsize=8)
    det=np.linspace(-8,8,1001)
    # Save locally resolved peak used in the figure separately.
    peak=[dict(detuning_over_hwhm=float(v),transmission=float(abs(local_transmission(.02,v*halfwidth(.02)))**2)) for v in det]
    write_csv('resolved_peak.csv',peak)
    ax[1].plot(det,[r['transmission'] for r in peak],color='#166b8a',label='Frozen matching response')
    ax[1].plot(det,1/(1+det*det),'--',color='#b16535',label='Local Lorentzian approximation')
    ax[1].axhline(.5,color='.7',lw=.6)
    ax[1].set(xlabel='Detuning / half-width',ylabel='Monochromatic transmission',title='Resolved resonance, omega r_e = 0.02')
    ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    fig.suptitle('Magnetic scalar signals: finite bandwidth\nPublished low-frequency response; frozen mouth coefficients across each packet',fontsize=12)
    save(fig,'01-bandwidth')
    fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    ax[0].plot([r['u_over_pi'] for r in pulses],[r['flux_over_first_energy'] for r in pulses],color='#166b8a')
    ax[0].set(xlabel='Retarded time / (pi ell)',ylabel='Energy flux / first-pulse energy (ell = 1)',title='Early delayed pulses\nomega r_e = 0.05, b = 10')
    for x in [.005,.02,.05]:
        rr=[r for r in collection if r['x']==x]
        ax[1].semilogx([r['arrival_over_pi_ell'] for r in rr],[r['fraction_of_eventual'] for r in rr],label=f'omega r_e = {x:g}')
    ax[1].axhline(.9,color='.6',ls=':',lw=.8)
    ax[1].set(xlabel='Last pulse-center time / (pi ell)',ylabel='Collected / eventual transmitted energy',title='Collection of the transmitted fraction\nNormalized by the small eventual signal')
    ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    fig.suptitle('Resumming repeated reflections gives a decaying pulse train\nConditional radial signal; ambient return paths excluded',fontsize=12)
    save(fig,'02-collection-time')
    fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    ax[0].loglog([r['x'] for r in delays],[r['sigma_t_over_ell_at_bandwidth_hwhm'] for r in delays],color='#166b8a',label='Packet sigma_t at b = half-width')
    ax[0].axhline(np.pi/eq['y'],color='#b16535',ls='--',label='Exterior d/ell at d/d0 = 0.01')
    ax[0].set(xlabel='Carrier omega r_e',ylabel='Packet time width / ell',title='Resonance resolution\nversus packet localization')
    ax[0].legend(fontsize=8)
    ax[1].loglog([r['q'] for r in physical],[r['omega_times_t90'] for r in physical],color='#166b8a',label='Orbital phase over 90% collection time')
    ax[1].axhline(1,color='.5',ls=':',label='One radian reference')
    ax[1].set(xlabel='Illustrative magnetic flux q (fixed g = 0.3)',ylabel='Omega_orbit t90',title='Exterior motion accumulated\nover the collection interval')
    ax[1].set_xlim(1e6,1e28)
    ax[1].set_xticks([1e6,1e10,1e14,1e18,1e22,1e26])
    ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    fig.suptitle('Additional experiment constraints, not a physical exclusion theorem\nRight: one quantum at 10% of binding energy; d/d0 = 0.01; other validity conditions pending',fontsize=12)
    save(fig,'03-experiment-constraints')
    sources=['scripts/magnetic_signals.py','scripts/run_magnetic_signals.py','scripts/magnetic_energy.py']
    summary=dict(author='Biswajit Jana',date='2026-11-25',model='Low-frequency neutral scalar matching; conditional radial energy; frozen mouth coefficients',
        band_cases=len(bands),independent_matching_cases=len(matching),band_audits=len(band_audit),time_audits=len(time_checks),
        compact_pulse_cases=len(compact),max_compact_pulse_difference=max(r['absolute_difference'] for r in compact),
        max_matching_difference=max(r['matching_error'] for r in matching),max_flux_error=max(r['flux_error'] for r in matching),
        max_band_fourier_difference=max(r['absolute_difference'] for r in band_audit),max_time_integral_difference=max(r['absolute_difference'] for r in time_checks),
        scenario=scenario[0],source_hashes={s:hashlib.sha256((ROOT/s).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for s in sources})
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
    paths=[*OUT.glob('*.csv'),OUT/'summary.json',*FIG.glob('*.png'),*FIG.glob('*.svg')]
    (OUT/'sha256.json').write_text(json.dumps({p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
