"""Reproduce gravity response checks; no course lab content. Runs in a temporary directory."""
import tempfile
import numpy as np
from fys2160_md import MDSimulation, Plummer

def diagnose(sim):
    dr=sim._x[:,None,:]-sim._x[None,:,:]
    r2=(dr*dr).sum(axis=2);np.fill_diagonal(r2,np.inf)
    specific=.5*(sim._v**2).sum(axis=1)-sim.G*(sim._mass[None,:]/np.sqrt(r2+sim.softening**2)).sum(axis=1)
    return int((specific>0).sum())

folder = tempfile.TemporaryDirectory(prefix='fys2160-gravity-validation-')
for label, dt, softening, removal in [('half-dt',.0025,.05,-.4),('less-softening',.0025,.025,-.4),('control',.005,.05,0)]:
    sim=MDSimulation(Plummer(N=256,G=1/256,softening=softening),pair_potential='gravity',G=1/256,softening=softening,ensemble='nve',timestep=dt,output_dir=folder.name)
    rows=[]
    for name,duration,rate in [('relax',20,0),('before',20,0),('remove',20,removal),('settle',20,0),('after',40,0)]:
        run=sim.run(steps=round(duration/dt),heat_rate=rate,thermo_every=round(.1/dt),trajectory_every=None,storage_name=label+'-'+name)
        d=run.thermo
        rows.append(d)
        if name in ('before','after'):
            print(label,name,'T',d.temperature.mean(),'E',d.total_energy.mean(),'rhalf',d.half_mass_radius.mean(),'virial',d.virial_ratio.mean(),'positive energy particles',diagnose(sim),flush=True)
    import pandas as pd
    d=pd.concat(rows);balance=d.total_energy-d.heat_added
    print(label,'energy balance range',balance.max()-balance.min(),flush=True)

folder.cleanup()
