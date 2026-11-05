"""SVPS Miniproject 1. Run: python solve.py. All quantities are dimensionless.

State: r, radial speed u, unwrapped angle theta, angular speed w,
       integrated thrust work W, integrated thrust torque J.
Initial position counts as crossing 1. Escape means the FIRST E=0.
"""
from pathlib import Path
import csv
import json
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent

def energy(z):
    r, u, _, w = z[:4]
    return (u*u + (r*w)**2)/2 - 1/r

def rhs(alpha, mode):
    if mode not in ('transverse', 'velocity'):
        raise ValueError(mode)
    def fun(t, z):
        r, u, theta, w = z[:4]
        vt = r*w
        speed = np.hypot(u, vt)
        if mode == 'transverse':
            ar, at = 0., alpha
        else:
            if speed == 0:
                raise ValueError('Velocity direction is undefined at zero speed')
            ar, at = alpha*u/speed, alpha*vt/speed
        return [u, r*w*w-1/r**2+ar, w, (at-2*u*w)/r,
                u*ar+vt*at, r*at]
    return fun

def escape(t, z):
    return energy(z)
escape.direction = 1
escape.terminal = True

def eighth(t, z):
    return z[2]-14*np.pi
eighth.direction = 1
eighth.terminal = False

def simulate(alpha=.005, mode='transverse', u0=0., w0=1.,
             rtol=1e-11, max_step=.2, method='DOP853'):
    z0 = [1., u0, 0., w0, 0., 0.]
    if energy(z0) >= 0:
        raise ValueError('Study requires an initially bound orbit')
    sol = solve_ivp(rhs(alpha, mode), (0, 10/alpha), z0,
                    method=method, rtol=rtol, atol=rtol*.01,
                    max_step=max_step, dense_output=True, events=(escape,eighth))
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError('Escape not resolved: '+sol.message)
    return sol

def summary(sol):
    t = float(sol.t_events[0][0])
    z = sol.y_events[0][0]
    return dict(escape_time=t, escape_radius=float(z[0]),
                escape_speed=float(np.hypot(z[1],z[0]*z[3])),
                escape_theta=float(z[2]), escape_x=float(z[0]*np.cos(z[2])),
                escape_y=float(z[0]*np.sin(z[2])),
                eighth_time=float(sol.t_events[1][0]) if len(sol.t_events[1]) else None)

def savecsv(name, rows):
    with (OUT/name).open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def savefig(fig, name):
    fig.tight_layout(pad=1.4)
    for ext in ('png','pdf','svg'):
        fig.savefig(OUT/f'{name}.{ext}', dpi=300)
    plt.close(fig)

def main():
    plt.rcParams.update({'font.family':['Arial','DejaVu Sans'], 'font.size':11,
        'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,
        'svg.fonttype':'none','axes.linewidth':1.1})
    base = simulate(); vel = simulate(mode='velocity')
    b, v = summary(base), summary(vel)
    t8 = b['eighth_time']; te = b['escape_time']
    z8 = base.sol(t8); ze = base.sol(te)
    ts = np.unique(np.r_[np.linspace(0,te,12001), t8,te])
    z = base.sol(ts); r,u,theta,w,W,J = z
    E = energy(z)
    checks = dict(max_work_residual=float(np.max(np.abs(E+.5-W))),
                  max_torque_residual=float(np.max(np.abs(r*r*w-1-J))),
                  escape_energy_residual=float(abs(energy(ze))),
                  eighth_angle_residual=float(abs(z8[2]-14*np.pi)))
    assert np.all(np.diff(E)>0)
    assert checks['max_work_residual'] < 1e-8
    assert checks['max_torque_residual'] < 1e-8
    crossings = []
    for k in range(8):
        t = 0. if k==0 else brentq(lambda t: base.sol(t)[2]-2*np.pi*k,0,te)
        q = base.sol(t)
        crossings.append(dict(crossing=k+1,theta_over_pi=2*k,time=t,radius=q[0],energy=energy(q)))
    savecsv('crossings.csv',crossings)
    savecsv('trajectory.csv',[dict(tau=t,r=q[0],u=q[1],theta=q[2],omega=q[3],
        x=q[0]*np.cos(q[2]),y=q[0]*np.sin(q[2]),energy=energy(q)) for t,q in zip(ts,z.T)])
    fig,ax = plt.subplots(figsize=(7,6))
    ax.plot(r*np.cos(theta),r*np.sin(theta),c='#0F4D92',lw=1.3,label='Transverse thrust')
    circle=np.linspace(0,2*np.pi,300)
    ax.plot(np.cos(circle),np.sin(circle),'--',c='0.65',lw=.8,label='Initial orbit (r = 1)')
    ax.plot(0,0,'k+',ms=8,label='Gravity centre')
    ax.scatter([z8[0]],[0],c='#B64342',s=40,zorder=4)
    ax.annotate(f'8th crossing\n'+r'$\tau$'+f' = {t8:.6f}',(z8[0],0),xytext=(20,20),textcoords='offset points',arrowprops={'arrowstyle':'->'})
    xe,ye=b['escape_x'],b['escape_y']
    ax.scatter([xe],[ye],marker='*',s=90,c='#B64342',zorder=4)
    ax.annotate(f'Escape (E = 0)\n'+r'$\tau$'+f' = {te:.6f}',(xe,ye),xytext=(-125,15),textcoords='offset points',arrowprops={'arrowstyle':'->'})
    ax.axhline(0,c='0.8',lw=.6); ax.axvline(0,c='0.8',lw=.6)
    ax.set(xlabel=r'$x=X/R_0$',ylabel=r'$y=Y/R_0$',title='Low-thrust orbit: required trajectory')
    ax.set_aspect('equal'); ax.margins(.18); ax.legend(loc='lower left',fontsize=9)
    savefig(fig,'figure1_trajectory')
    fig,ax=plt.subplots(figsize=(7,4.3))
    ax.plot(ts,E,c='#0F4D92',lw=2)
    ax.axhline(0,c='0.4',ls='--',lw=1)
    ax.scatter([t8,te],[energy(z8),0],c='#B64342',zorder=3)
    ax.annotate('8th crossing',(t8,energy(z8)),xytext=(8,-24),textcoords='offset points')
    ax.annotate(f'Escape: {te:.6f}',(te,0),xytext=(-120,-55),textcoords='offset points',arrowprops={'arrowstyle':'->','color':'0.4'})
    ax.set(xlabel=r'Dimensionless time $\tau=t/T_0$',ylabel=r'Specific energy $\mathcal{E}=\epsilon/(GM/R_0)$',
           title='Energy increases from -0.5 to 0')
    savefig(fig,'figure2_energy')
    convergence=[]
    for method,tol,step in [('DOP853',1e-8,1.),('DOP853',1e-10,.5),('DOP853',1e-12,.1),('RK45',1e-11,.1)]:
        s=summary(simulate(rtol=tol,max_step=step,method=method))
        convergence.append(dict(method=method,rtol=tol,max_step=step,t8=s['eighth_time'],tesc=s['escape_time'],
                                delta_t8=s['eighth_time']-t8,delta_tesc=s['escape_time']-te))
    savecsv('convergence.csv',convergence)
    # Independent Cartesian formulation, including continuous transverse thrust.
    def cart(t,q):
        x,y,vx,vy=q; rr=np.hypot(x,y)
        return [vx,vy,-x/rr**3-.005*y/rr,-y/rr**3+.005*x/rr]
    def ce(t,q):
        return (q[2]**2+q[3]**2)/2-1/np.hypot(q[0],q[1])
    ce.direction=1; ce.terminal=True
    cs=solve_ivp(cart,(0,300),[1,0,0,1],events=ce,dense_output=True,
                 method='DOP853',rtol=1e-12,atol=1e-14,max_step=.1)
    ct8=brentq(lambda t:cs.sol(t)[1],t8-.1,t8+.1)
    assert cs.sol(ct8)[0]>0
    checks.update(cartesian_escape_difference=float(cs.t_events[0][0]-te),cartesian_t8_difference=float(ct8-t8))
    assert abs(checks['cartesian_escape_difference'])<1e-6
    # Zero-thrust analytic circular orbit regression.
    zero=solve_ivp(rhs(0,'transverse'),(0,14*np.pi),[1,0,0,1,0,0],rtol=1e-12,atol=1e-14,method='DOP853')
    checks['zero_thrust_state_error']=float(np.max(np.abs(zero.y[:4,-1]-[1,0,14*np.pi,1])))
    assert checks['zero_thrust_state_error']<1e-8
    sweep=[]
    initial=[('circular',0,1),('outward',.2,1),('inward',-.2,1),('slow',0,.9),('fast',0,1.1)]
    for label,u0,w0 in initial:
        for alpha in [.001,.0025,.005,.01,.02]:
            a=summary(simulate(alpha,'transverse',u0,w0,max_step=.5))
            c=summary(simulate(alpha,'velocity',u0,w0,max_step=.5))
            err=abs(a['escape_time']-c['escape_time'])/c['escape_time']
            sweep.append(dict(initial=label,u0=u0,w0=w0,alpha=alpha,
                transverse_escape=a['escape_time'],velocity_escape=c['escape_time'],
                relative_error_percent=err*100,pass_1percent=bool(err<=.01)))
    savecsv('extension_sweep.csv',sweep)
    fig,axes=plt.subplots(1,2,figsize=(10,4.2))
    for sol,label,c in [(base,'Transverse','#0F4D92'),(vel,'Velocity aligned','#B64342')]:
        tt=np.linspace(0,sol.t[-1],3000)
        axes[0].plot(tt,energy(sol.sol(tt)),label=label,c=c,lw=2)
    axes[0].set(xlabel=r'$\tau$',ylabel=r'$\mathcal{E}$',title='A: energy history'); axes[0].legend()
    for label,_,_ in initial:
        rows=[s for s in sweep if s['initial']==label]
        axes[1].plot([s['alpha'] for s in rows],[s['relative_error_percent'] for s in rows],marker='o',ms=4,label=label)
    axes[1].axhline(1,c='k',ls='--',lw=1,label='1% tolerance')
    axes[1].set(xscale='log',xlabel=r'Thrust parameter $\alpha$',ylabel='Escape-time error (%)',title='A: sensitivity to initial state')
    axes[1].legend(fontsize=8,ncol=2)
    savefig(fig,'figure3_extension')
    results=dict(baseline=b,velocity_aligned=v,eighth_radius=float(z8[0]),eighth_energy=float(energy(z8)),
        checks=checks,versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__))
    (OUT/'results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results,indent=2))

if __name__=='__main__':
    main()
