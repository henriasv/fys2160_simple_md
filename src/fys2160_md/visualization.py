"""Bounded notebook visualization; no external JavaScript or continuous render loop."""
import json
import time
import uuid
import numpy as np


class LiveView:
    """Throttled x–y projection while C advances the simulation in batches."""
    def __init__(self):
        from IPython import get_ipython
        from IPython.display import display, HTML
        if get_ipython() is None:
            raise ValueError('show=True requires IPython/Jupyter; use Run.view() for saved trajectories.')
        self.handle = display(HTML('<p>Starting molecular dynamics…</p>'), display_id=True)
        self.last = -np.inf
        self.identifier = "md-live-" + uuid.uuid4().hex

    def update(self, sim, *, force=False):
        now = time.monotonic()
        if not force and now-self.last < .25:
            return
        from IPython.display import HTML
        self.last = now
        # At most 1200 circles; recording is independent of this display limit.
        positions = sim._x[:1200]/sim._box
        colors = ('blue', 'orange')
        gradients = ''.join(
            f'<radialGradient id="{self.identifier}-{name}" cx="30%" cy="25%" r="78%">'
            f'<stop offset="0" stop-color="{light}"/><stop offset=".55" stop-color="{middle}"/>'
            f'<stop offset="1" stop-color="{dark}"/></radialGradient>'
            for name, light, middle, dark in [('blue','#e7f6ff','#3698bc','#153d56'),
                                              ('orange','#fff0d8','#de9148','#754019')])
        circles = ''.join(
            f'<circle cx="{20+480*x:.2f}" cy="{500-480*y:.2f}" r="3.5" '
            f'fill="url(#{self.identifier}-{colors[i%2 if sim.model != "atomic" else 0]})"/>'
            for i,(x,y,z) in enumerate(positions))
        obs = sim.observables
        caption = f'Step {sim.step} · T* = {obs["temperature"]:.3f} · P* = {obs["pressure"]:.5f} · x–y projection'
        if len(sim._x)>1200:
            caption += ' · showing first 1200 atoms'
        self.handle.update(HTML(f'<div><svg viewBox="0 0 520 520" width="420" role="img" aria-label="Live atom positions"><defs>{gradients}</defs><rect x="20" y="20" width="480" height="480" fill="#f5f8fa" stroke="#526879"/>{circles}</svg><p>{caption}</p></div>'))


def trajectory_player(run, *, max_frames=150, max_atoms=1500, projection='orthographic', zoom=1.0):
    from IPython.display import HTML
    zoom = float(zoom)
    if not np.isfinite(zoom) or not .25 <= zoom <= 3:
        raise ValueError('zoom must be between 0.25 and 3.')
    if projection not in ('orthographic', 'perspective'):
        raise ValueError('projection must be orthographic or perspective.')
    if int(max_frames) != max_frames or max_frames < 2 or int(max_atoms) != max_atoms or max_atoms < 2:
        raise ValueError('max_frames and max_atoms must be integers >= 2.')
    paths = sorted((run.path/'frames').glob('*.npz'))
    if not paths:
        raise ValueError('No trajectory frames were saved.')
    selected = np.unique(np.linspace(0, len(paths)-1, min(int(max_frames), len(paths))).astype(int))
    frames = []
    for index in selected:
        with np.load(paths[index], allow_pickle=False) as f:
            frames.append(dict(x=np.round(f['positions'][:int(max_atoms)], 5).tolist(),
                               box=f['box'].tolist(), time=float(f['time']), step=int(f['step'])))
    data = json.dumps(dict(frames=frames, projection=projection, zoom=zoom, molecular=run.metadata['configuration']['model'] != 'atomic'))
    identifier = 'md-'+uuid.uuid4().hex
    note = f'{len(frames)} of {len(paths)} saved frames; {len(frames[0]["x"])} of {run.metadata["atoms"]} atoms. Drag to rotate; scroll over the view to zoom. Arrow keys also rotate the focused view. Display sampling does not change saved data.'
    template = '''<div id="__ID__" style="max-width:720px;font:14px system-ui;color:#233744">
<canvas tabindex="0" width="720" height="600" style="width:100%;background:#f4f7fa;touch-action:none" aria-label="Molecular trajectory, drag to rotate"></canvas>
<div style="display:flex;gap:12px;align-items:center;margin:10px 0"><button type="button">Play</button><button type="button" data-reset>Reset view</button><input aria-label="Trajectory frame" type="range" min="0" value="0" style="flex:1"><output></output></div><div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap"><label>Projection <select aria-label="Projection"><option value="orthographic">Orthographic</option><option value="perspective">Perspective</option></select></label><label style="display:flex;align-items:center;gap:8px">Zoom <button type="button" aria-label="Zoom out" data-zoom-out>−</button><input type="range" aria-label="Zoom" min="25" max="300" step="1" value="100" style="width:130px"><button type="button" aria-label="Zoom in" data-zoom-in>+</button><span data-zoom-value>100%</span></label></div><p>__NOTE__</p>
<script>(function(){
const root=document.getElementById('__ID__'),data=__DATA__,frames=data.frames;
const canvas=root.querySelector('canvas'),ctx=canvas.getContext('2d'),button=root.querySelector('button'),slider=root.querySelector('input'),label=root.querySelector('output');
const projection=root.querySelector('select');projection.value=data.projection;
const zoomSlider=root.querySelector('[aria-label="Zoom"]'),zoomLabel=root.querySelector('[data-zoom-value]');
let zoom=data.zoom;zoomSlider.value=Math.round(100*zoom);zoomLabel.textContent=Math.round(100*zoom)+'%';
slider.max=frames.length-1;const extent=Math.max(...frames.flatMap(f=>f.box));
function normalize(p,box){return p.map((v,d)=>(v-.5*box[d])/extent+.5);}
let drag=null,timer=null;
const cy=Math.cos(.35),sy=Math.sin(.35),cp=Math.cos(.3),sp=Math.sin(.3);
const initial=[[cy,0,-sy],[-sy*sp,cp,-cy*sp],[sy*cp,sp,cy*cp]];
let rotation=initial.map(row=>row.slice());
// Apply camera-axis rotations, so the near face follows the pointer at any orientation.
function rotateView(dx,dy){
  const cY=Math.cos(dx),sY=Math.sin(dx),cX=Math.cos(dy),sX=Math.sin(dy);
  for(let k=0;k<3;k++){
    const x=rotation[0][k],y=rotation[1][k],z=rotation[2][k];
    const xx=cY*x+sY*z,zz=-sY*x+cY*z;
    rotation[0][k]=xx;rotation[1][k]=cX*y-sX*zz;rotation[2][k]=sX*y+cX*zz;
  }
}
function project(p){const centered=p.map(v=>v-.5);
  const r=rotation.map(row=>row.reduce((sum,v,i)=>sum+v*centered[i],0));
  const factor=projection.value==='perspective'?2.5/(2.5-r[2]):1;
  return [360+300*zoom*factor*r[0],300-300*zoom*factor*r[1],r[2],factor];
}
// Build two shaded billboard textures once; drawing particles stays inexpensive.
function sphereTexture(stops){
  const sprite=document.createElement('canvas');sprite.width=sprite.height=64;
  const c=sprite.getContext('2d'),g=c.createRadialGradient(22,20,1,30,30,33);
  stops.forEach(([offset,color])=>g.addColorStop(offset,color));
  c.fillStyle=g;c.beginPath();c.arc(32,32,30,0,2*Math.PI);c.fill();return sprite;
}
const sprites=[sphereTexture([[0,'#f0fbff'],[.3,'#9cd7ed'],[.65,'#3494ba'],[1,'#153d56']]),
               sphereTexture([[0,'#fff5df'],[.3,'#f6c186'],[.65,'#d78238'],[1,'#754019']])];
function line(a,b,color,width){ctx.strokeStyle=color;ctx.lineWidth=width;ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();}
function draw(){const f=frames[+slider.value];ctx.clearRect(0,0,720,600);
for(let a=0;a<8;a++)for(let d=0;d<3;d++){const b=a^(1<<d);if(a<b)line(project(normalize([a&1,(a>>1)&1,(a>>2)&1].map((v,d)=>v*f.box[d]),f.box)),project(normalize([b&1,(b>>1)&1,(b>>2)&1].map((v,d)=>v*f.box[d]),f.box)),'#a8b8c4',1);}
const positions=f.x.map(p=>p.map((v,d)=>v/f.box[d]));const points=f.x.map(p=>project(normalize(p,f.box)));
if(data.molecular)for(let i=0;i+1<points.length;i+=2){if(positions[i].every((v,d)=>Math.abs(v-positions[i+1][d])<.5))line(points[i],points[i+1],'#7d8b96',2);}
points.map((p,i)=>({p,i})).sort((a,b)=>a.p[2]-b.p[2]).forEach(({p,i})=>{const radius=4.5*zoom*p[3];ctx.drawImage(sprites[data.molecular&&i%2?1:0],p[0]-radius,p[1]-radius,2*radius,2*radius);});
label.textContent='Step '+f.step+' · t* '+f.time.toFixed(2);}
function stop(){clearInterval(timer);timer=null;button.textContent='Play';}
button.onclick=()=>{if(timer){stop();return;}button.textContent='Pause';timer=setInterval(()=>{if(!root.isConnected){stop();return;}slider.value=(+slider.value+1)%frames.length;draw();},100);};
projection.onchange=draw;slider.oninput=draw;canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);};
canvas.onpointermove=e=>{if(!drag)return;rotateView((e.clientX-drag[0])*.01,(e.clientY-drag[1])*.01);drag=[e.clientX,e.clientY];draw();};canvas.onpointerup=canvas.onpointercancel=()=>drag=null;
function setZoom(value){zoom=Math.max(.25,Math.min(3,value));zoomSlider.value=Math.round(100*zoom);zoomLabel.textContent=Math.round(100*zoom)+'%';draw();}
zoomSlider.oninput=()=>setZoom(+zoomSlider.value/100);
root.querySelector('[data-zoom-out]').onclick=()=>setZoom(zoom/1.2);
root.querySelector('[data-zoom-in]').onclick=()=>setZoom(zoom*1.2);
canvas.addEventListener('wheel',e=>{e.preventDefault();setZoom(zoom*Math.exp(-e.deltaY*.001));},{passive:false});
root.querySelector('[data-reset]').onclick=()=>{rotation=initial.map(row=>row.slice());setZoom(data.zoom);};
canvas.onkeydown=e=>{const delta={ArrowRight:[.08,0],ArrowLeft:[-.08,0],ArrowDown:[0,.08],ArrowUp:[0,-.08]}[e.key];if(delta){e.preventDefault();rotateView(...delta);draw();}};
draw();})();</script></div>'''
    return HTML(template.replace('__ID__',identifier).replace('__NOTE__',note).replace('__DATA__',data))
