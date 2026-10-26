/* Shared 3D renderer for live widgets and self-contained saved players.
   OMOVI (github.com/andeplane/omovi, by Anders Hafreager) renders the scene when it
   can be loaded from jsDelivr; the built-in WebGL renderer draws first and remains
   the fallback offline. */
const OMOVI_URL='https://cdn.jsdelivr.net/npm/omovi@0.25.0/+esm';
const THREE_URL='https://cdn.jsdelivr.net/npm/three@0.181.2/+esm';  // the three.js that omovi imports
const PALETTE=[[41,145,184],[214,122,48],[64,166,102],[153,102,191],[204,77,89],[128,140,153],[230,232,235]];
// omovi 0.25 caches its particle/bond materials per module, so visualizers in one
// page would share colours and radii. A distinct URL per view gives each its own
// module instance (three.js stays shared).
function loadOmovi(){
  const view=globalThis.__fys2160OmoviViews=(globalThis.__fys2160OmoviViews||0)+1;
  return Promise.all([import(OMOVI_URL+'?view='+view),import(THREE_URL)]);
}
// The colour behind an element: the first non-transparent ancestor background, so
// views blend into light and dark pages alike.
function pageBackground(node){
  for(let e=node;e&&e.nodeType===1;e=e.parentElement){
    const m=getComputedStyle(e).backgroundColor.match(/rgba?\(([^)]+)\)/);
    if(m){const [r,g,b,a=1]=m[1].split(',').map(Number);if(a>0)return [r,g,b];}
  }
  return window.matchMedia?.('(prefers-color-scheme: dark)').matches?[30,30,30]:[255,255,255];
}
const ICONS={
  play:'<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M8 5.2v13.6L18.6 12z"/></svg>',
  pause:'<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M6.5 5h4v14h-4zm7 0h4v14h-4z"/></svg>',
  reset:'<svg viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"><path d="M4.5 12a7.5 7.5 0 1 0 2.2-5.3"/><path d="M4.5 4.5v4h4"/></svg>'
};
function escapeHTML(text){return String(text).replace(/[&<>"']/g,c=>'&#'+c.charCodeAt(0)+';');}
// The controls take their colours from the surrounding text colour, so they suit
// light and dark pages; the accent is the first palette colour.
const VIEW_STYLE=`
[data-md-view]{--md-accent:#2991b8;--md-line:color-mix(in srgb,currentColor 16%,transparent);--md-soft:color-mix(in srgb,currentColor 5%,transparent);--md-muted:color-mix(in srgb,currentColor 65%,transparent)}
[data-md-view] [hidden]{display:none!important}
[data-md-view] .md-card{border:1px solid var(--md-line)}
[data-md-view] .md-card>[tabindex]:focus-visible{outline:2px solid var(--md-accent);outline-offset:-2px}
[data-md-view] .md-bar{border-top:1px solid var(--md-line);background:var(--md-soft);padding:6px 8px}
[data-md-view] .md-row{display:flex;align-items:center;gap:4px;flex-wrap:wrap}
[data-md-view] .md-btn,[data-md-view] select{all:unset;box-sizing:border-box;height:30px;border-radius:7px;cursor:pointer;color:inherit;font:600 13px/30px system-ui,sans-serif}
[data-md-view] .md-btn{display:inline-flex;align-items:center;justify-content:center;min-width:30px;padding:0 7px}
[data-md-view] .md-btn:hover,[data-md-view] select:hover{background:var(--md-line)}
[data-md-view] .md-btn:focus-visible,[data-md-view] select:focus-visible,[data-md-view] input:focus-visible{outline:2px solid var(--md-accent);outline-offset:1px}
[data-md-view] .md-btn svg{width:16px;height:16px;fill:currentColor}
[data-md-view] .md-btn svg[stroke]{fill:none}
[data-md-view] .md-play{width:32px;height:32px;padding:0;border-radius:50%;background:var(--md-accent);color:#fff}
[data-md-view] .md-play:hover{background:color-mix(in srgb,var(--md-accent) 82%,#000)}
[data-md-view] input[type=range]{flex:1 1 160px;min-width:40px;margin:0 6px;accent-color:var(--md-accent);cursor:pointer}
[data-md-view] .md-count,[data-md-view] output{font:12px/1.4 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-variant-numeric:tabular-nums;color:var(--md-muted)}
[data-md-view] .md-count{min-width:9ch;text-align:right;margin-right:4px}
[data-md-view] .md-tools{display:flex;align-items:center;gap:2px;margin-left:auto}
[data-md-view] select{padding:0 24px 0 9px;font-weight:500;background:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 10 6'%3E%3Cpath d='M1 1l4 4 4-4' fill='none' stroke='%23888' stroke-width='1.6'/%3E%3C/svg%3E") no-repeat right 8px center/10px}
[data-md-view] select option{color:#222;background:#fff}
[data-md-view] .md-info{gap:4px 12px;padding-left:2px}
[data-md-view] [data-playback]{flex-wrap:nowrap}
[data-md-view] [data-playback]+.md-info{margin-top:2px}
[data-md-view] .md-legend{display:flex;flex-wrap:wrap;gap:10px;font-size:12px;color:var(--md-muted)}
[data-md-view] .md-chip{display:inline-flex;align-items:center;gap:5px}
[data-md-view] .md-chip::before{content:"";width:10px;height:10px;border-radius:50%;background:var(--c);box-shadow:inset 0 0 0 1px rgba(0,0,0,.25)}
[data-md-view] .md-status{display:inline-flex;align-items:center;gap:7px;height:30px;padding:0 11px 0 9px;border-radius:15px;background:var(--md-line);font:600 12px system-ui,sans-serif}
[data-md-view] .md-dot{width:8px;height:8px;border-radius:50%;background:#8a949c}
[data-md-view] [data-state=running] .md-dot{background:#3aa76d;animation:md-pulse 1.6s ease-in-out infinite}
[data-md-view] [data-state=interrupted] .md-dot{background:#d99a2b}
[data-md-view] [data-state=failed] .md-dot{background:#d0504a}
@keyframes md-pulse{50%{opacity:.3}}
@media (prefers-reduced-motion:reduce){[data-md-view] .md-dot{animation:none}}
[data-md-view] [data-note]{margin:8px 2px 0;font-size:12.5px;line-height:1.5;color:var(--md-muted)}
`;
const VIEW_MARKUP=`<style>${VIEW_STYLE}</style>
<div class="md-card">
<canvas tabindex="0" width="720" height="600" style="display:block;width:100%;height:auto;aspect-ratio:6/5;background:transparent;touch-action:none" aria-label="Molecular trajectory, drag to rotate"></canvas>
<div class="md-bar">
<div class="md-row" data-playback>
<button type="button" class="md-btn md-play" data-play aria-label="Play" title="Play (space)">${ICONS.play}</button>
<input type="range" aria-label="Trajectory frame" min="0" value="0">
<span class="md-count" data-count></span>
</div>
<div class="md-row md-info">
<span class="md-status" data-status hidden><span class="md-dot"></span><span data-status-text></span></span>
<output></output><span class="md-legend" data-legend></span>
<span class="md-tools">
<select data-speed aria-label="Playback speed" title="Playback speed (&lt; and &gt; keys)"><option value=".25">0.25\u00d7</option><option value=".5">0.5\u00d7</option><option value="1" selected>1\u00d7</option><option value="2">2\u00d7</option><option value="4">4\u00d7</option><option value="8">8\u00d7</option></select>
<button type="button" class="md-btn" data-zoom-out aria-label="Zoom out" title="Zoom out">−</button>
<button type="button" class="md-btn" data-zoom-in aria-label="Zoom in" title="Zoom in">+</button>
<button type="button" class="md-btn" data-reset aria-label="Reset view" title="Reset view">${ICONS.reset}</button>
<select data-projection aria-label="Projection" title="Projection"><option value="perspective">Perspective</option><option value="orthographic">Orthographic</option></select>
<button type="button" class="md-btn" data-help aria-label="Controls and details" aria-expanded="false" title="Controls and details">?</button>
</span>
</div>
</div>
</div>
<p data-note hidden></p>`;
function createMDView(el,data){
// Own the sizing element: widget hosts may manage/reset styles on their `el`.
// A fixed CSS width scales with editor zoom; max-width only shrinks narrow panes.
const root=document.createElement('div'),frames=data.frames;
root.dataset.mdView='';
root.style.cssText='display:block;box-sizing:border-box;width:720px;max-width:100%;min-width:0;font:14px system-ui;color:inherit';
el.replaceChildren(root);
root.innerHTML=VIEW_MARKUP;
const speciesNames=[...new Set(data.species||[])];
const speciesColors=['blue','orange','green','purple','red','grey','white'];
// Conventional element colours for water; other species use the palette in order.
const elementKinds={O:4,H:6};
const kindOf=name=>name in elementKinds?elementKinds[name]:speciesNames.indexOf(name)%6;
const note=root.querySelector('[data-note]'),help=root.querySelector('[data-help]');
note.textContent=data.note+(speciesNames.length>1?' Species: '+speciesNames.map(name=>name+' ('+speciesColors[kindOf(name)]+')').join(', ')+'.':'');
help.onclick=()=>{note.hidden=!note.hidden;help.setAttribute('aria-expanded',String(!note.hidden));};
// Species colour chips next to the readout; omovi may replace the colours.
const legend=root.querySelector('[data-legend]');
function showLegend(colorOf){
  legend.innerHTML=speciesNames.length>1?speciesNames.map(name=>'<span class="md-chip" style="--c:rgb('+colorOf(name).join(',')+')">'+escapeHTML(name)+'</span>').join(''):'';
}
showLegend(name=>PALETTE[kindOf(name)]);
let pendingDraw=0;
const canvas=root.querySelector('canvas'),renderer=createMDSphereRenderer(canvas),button=root.querySelector('[data-play]'),slider=root.querySelector('input'),label=root.querySelector('output');
const count=root.querySelector('[data-count]'),speedSelect=root.querySelector('[data-speed]');
const status=root.querySelector('[data-status]'),statusText=root.querySelector('[data-status-text]');
const projection=root.querySelector('[data-projection]');projection.value=data.projection;
let omovi=null,disposed=false;
let zoom=data.zoom;
if(data.live){
  root.querySelector('[data-playback]').hidden=true;
  button.hidden=slider.hidden=count.hidden=speedSelect.hidden=true;status.hidden=false;
  canvas.setAttribute('aria-label','Live molecular simulation, drag to rotate');
}
slider.max=frames.length-1;const extent=Math.max(...frames.flatMap(f=>f.box));
function normalize(p,box){return p.map((v,d)=>(v-.5*box[d])/extent+.5);}
const pointers=new Map();
let pan=[0,0],timer=null;
const cy=Math.cos(.35),sy=Math.sin(.35),cp=Math.cos(.3),sp=Math.sin(.3);
const initial=[[cy,0,-sy],[-sy*sp,cp,-cy*sp],[sy*cp,sp,cy*cp]];
// Orthographic views open face-on to the xz plane with z up, or with z across
// for a box elongated along z (a slab). Rows: screen right, screen up, towards viewer.
function frontView(){
  const b=frames[+slider.value].box;
  return b[2]===Math.max(...b)&&b[2]>1.5*Math.min(...b)?[[0,0,1],[1,0,0],[0,1,0]]:[[1,0,0],[0,0,1],[0,-1,0]];
}
function homeView(){return (projection.value==='orthographic'?frontView():initial).map(row=>row.slice());}
let rotation=homeView();
// Apply camera-axis rotations, so the near face follows the pointer at any orientation.
function rotateView(dx,dy){
  const cY=Math.cos(dx),sY=Math.sin(dx),cX=Math.cos(dy),sX=Math.sin(dy);
  for(let k=0;k<3;k++){
    const x=rotation[0][k],y=rotation[1][k],z=rotation[2][k];
    const xx=cY*x+sY*z,zz=-sY*x+cY*z;
    rotation[0][k]=xx;rotation[1][k]=cX*y-sX*zz;rotation[2][k]=sX*y+cX*zz;
  }
}
function cameraPosition(p){
  const centered=p.map(v=>v-.5);
  return rotation.map(row=>row.reduce((sum,v,i)=>sum+v*centered[i],0));
}
// Keep projection, particles and bonds in the same logical coordinate system.
// Rebuild the backing bitmap at the actual display resolution after editor zoom
// or pane resizing; changing its attributes must not change the CSS aspect ratio.
function draw(){
if(omovi){omovi.show(frames[+slider.value]);updateLabel(frames[+slider.value]);return;}
const bounds=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1;
const width=Math.max(1,Math.round(bounds.width*dpr));
const height=Math.max(1,Math.round(bounds.height*dpr));
if(canvas.width!==width||canvas.height!==height){canvas.width=width;canvas.height=height;}
const f=frames[+slider.value],edges=[],bonds=[];
const corners=Array.from({length:8},(_,a)=>cameraPosition(normalize([a&1,(a>>1)&1,(a>>2)&1].map((v,d)=>v*f.box[d]),f.box)));
for(let a=0;a<8;a++)for(let d=0;d<3;d++){const b=a^(1<<d);if(a<b&&data.boundary!=="open")edges.push(corners[a],corners[b]);}
const positions=f.x.map(p=>cameraPosition(normalize(p,f.box)));
const pairs=data.bonds||(data.molecular?Array.from({length:Math.floor(positions.length/2)},(_,k)=>[2*k,2*k+1]):[]);
for(const [i,j] of pairs){
  if(j<positions.length&&f.x[i].every((v,d)=>Math.abs(v-f.x[j][d])<.5*f.box[d]))bonds.push(positions[i],positions[j]);
}
renderer.draw({positions,edges,bonds,radius:.5*(typeof data.sigma==='number'?data.sigma:1)/extent,
radii:positions.map((_,i)=>.5*(typeof data.sigma==='object'?data.sigma[data.species[i]]:(data.sigma??1))/extent),
kinds:positions.map((_,i)=>speciesNames.length>0?kindOf(data.species[i]):(data.molecular?i%2:0)),zoom,pan,perspective:projection.value==='perspective',molecular:data.molecular});
updateLabel(f);}
function updateLabel(f){
  const u=data.units;
  if(!data.live)count.textContent=(+slider.value+1)+' / '+frames.length;
  label.textContent='Step '+f.step+(u?' · t '+(f.time*u.time_scale).toFixed(2)+' '+u.time:' · t* '+f.time.toFixed(2))
    +(data.live?(u?' · T '+f.temperature.toFixed(1)+' '+u.temperature:' · T* '+f.temperature.toFixed(3))
      +(Number.isFinite(f.pressure)?(u?' · P '+f.pressure.toPrecision(4)+' '+u.pressure:' · P* '+f.pressure.toPrecision(4)):''):'');
}
function showPlaying(playing){
  const name=playing?'Pause':'Play';
  button.innerHTML=playing?ICONS.pause:ICONS.play;button.setAttribute('aria-label',name);button.title=name+' (space)';
}
function stop(){clearInterval(timer);timer=null;showPlaying(false);}
// 10 frames per second at 1x.
let speed=1;
function start(){
  clearInterval(timer);showPlaying(true);
  timer=setInterval(()=>{if(!root.isConnected){stop();return;}slider.value=(+slider.value+1)%frames.length;draw();},100/speed);
}
function togglePlay(){if(timer)stop();else start();}
button.onclick=togglePlay;
speedSelect.onchange=()=>{speed=+speedSelect.value;if(timer)start();};
// < and > step the playback speed down and up, as in video players.
function stepSpeed(step){
  const values=[...speedSelect.options].map(o=>o.value),i=values.indexOf(speedSelect.value)+step;
  if(i>=0&&i<values.length){speedSelect.value=values[i];speedSelect.onchange();}
}
function setStatus(value){
  statusText.textContent=value.charAt(0).toUpperCase()+value.slice(1);status.dataset.state=value;
}
if(data.live)setStatus('running');
// In a focused saved view, space plays and pauses and < > change the speed; live
// views have no playback.
function playKey(e){
  if(data.live)return false;
  if(e.key==='<'||e.key==='>'){e.preventDefault();stepSpeed(e.key==='>'?1:-1);return true;}
  if(e.key!==' ')return false;
  e.preventDefault();if(!e.repeat)togglePlay();
  return true;
}
projection.onchange=()=>{if(omovi)omovi.setProjection(projection.value);draw();};slider.oninput=draw;
slider.onkeydown=playKey;
canvas.onpointerdown=e=>{
  if(e.button!==0&&e.button!==2)return;
  e.preventDefault();canvas.setPointerCapture(e.pointerId);
  pointers.set(e.pointerId,{x:e.clientX,y:e.clientY,pan:e.button===2||e.shiftKey});
};
function centroid(){
  const points=[...pointers.values()];
  return [points.reduce((s,p)=>s+p.x,0)/points.length,points.reduce((s,p)=>s+p.y,0)/points.length];
}
canvas.onpointermove=e=>{
  const previous=pointers.get(e.pointerId);if(!previous)return;
  const before=centroid(),dx=e.clientX-previous.x,dy=e.clientY-previous.y;
  previous.x=e.clientX;previous.y=e.clientY;
  const bounds=canvas.getBoundingClientRect();
  if(pointers.size>1||previous.pan){
    const after=centroid();
    pan[0]+=(after[0]-before[0])*720/bounds.width;
    pan[1]+=(after[1]-before[1])*600/bounds.height;
  }else rotateView(dx*720/bounds.width*.01,dy*600/bounds.height*.01);
  draw();
};
canvas.onpointerup=canvas.onpointercancel=canvas.onlostpointercapture=e=>pointers.delete(e.pointerId);
canvas.oncontextmenu=e=>e.preventDefault();
function setZoom(value){zoom=Math.max(.25,Math.min(3,value));draw();}
root.querySelector('[data-zoom-out]').onclick=()=>setZoom(zoom/1.2);
root.querySelector('[data-zoom-in]').onclick=()=>setZoom(zoom*1.2);
canvas.addEventListener('wheel',e=>{e.preventDefault();setZoom(zoom*Math.exp(e.deltaY*.001));},{passive:false});
root.querySelector('[data-reset]').onclick=()=>{rotation=homeView();pan=[0,0];if(omovi)omovi.resetCamera();setZoom(data.zoom);};
canvas.onkeydown=e=>{if(playKey(e))return;const delta={ArrowRight:[.08,0],ArrowLeft:[-.08,0],ArrowDown:[0,.08],ArrowUp:[0,-.08]}[e.key];if(delta){e.preventDefault();rotateView(...delta);draw();}};
function requestDraw(){
  if(!pendingDraw)pendingDraw=requestAnimationFrame(()=>{pendingDraw=0;draw();});
}
const resizeObserver=new ResizeObserver(requestDraw);
resizeObserver.observe(canvas);
// A display/zoom change can change DPR without changing the CSS dimensions.
let resolution;
function watchResolution(){
  if(resolution)resolution.removeEventListener('change',watchResolution);
  resolution=window.matchMedia(`(resolution: ${window.devicePixelRatio||1}dppx)`);
  resolution.addEventListener('change',watchResolution);
  requestDraw();
}
watchResolution();
draw();
// Per-atom omovi types: one per species (or per alternating diatomic atom), with
// the same colours and display radii as the built-in renderer.
function atomTypes(){
  const count=frames[0].x.length,sigma=i=>typeof data.sigma==='object'?data.sigma[data.species[i]]:(data.sigma??1);
  const types=[],styles=new Map();
  for(let i=0;i<count;i++){
    const key=speciesNames.length>0?speciesNames.indexOf(data.species[i]):(data.molecular?i%2:0);
    const kind=speciesNames.length>0?kindOf(data.species[i]):key;
    if(!styles.has(key+1))styles.set(key+1,{color:PALETTE[kind%7],radius:.5*sigma(i)});
    types.push(key+1);
  }
  return {types,styles};
}
function startOmovi([lib,THREE]){
  if(disposed)return;
  const container=document.createElement('div');
  container.style.cssText='display:block;width:100%;aspect-ratio:6/5;touch-action:none;position:relative;overflow:hidden';
  container.setAttribute('aria-label',canvas.getAttribute('aria-label'));
  // Focusable, and focused by a click or touch, so that space toggles playback.
  container.tabIndex=0;
  container.addEventListener('pointerdown',()=>container.focus({preventScroll:true}),true);
  container.addEventListener('keydown',playKey);
  canvas.replaceWith(container);  // in the page first: omovi measures its element
  const visualizer=new lib.Visualizer({domElement:container});
  container.omovi=visualizer;container.omoviLibrary=lib;  // for inspection from the browser console
  // omovi gives its canvas a 640x480 CSS size within percentage limits, and measures
  // that canvas when the element has no size yet; fill the view instead.
  visualizer.canvas.style.cssText='display:block;width:100%;height:100%';
  renderer.dispose();
  const lost=canvas.getContext('webgl2')?.getExtension('WEBGL_lose_context');if(lost)lost.loseContext();
  const {types,styles}=atomTypes();
  // omovi looks colours and radii up by particle id (not by type): set them per atom.
  let smallest=Infinity,largest=0,mean=0;
  // In real units, species named after elements (O, H, ...) use omovi's element
  // colours, as in Atomify; reduced-unit labels such as "B" are not boron.
  const element=i=>data.units&&speciesNames.length>0?lib.AtomTypes?.[data.species[i]]:null;
  types.forEach((t,i)=>{
    const style=styles.get(t),own=element(i)?.color;
    visualizer.setColor(i,own||{r:style.color[0],g:style.color[1],b:style.color[2]});
    visualizer.setRadius(i,style.radius);
    smallest=Math.min(smallest,style.radius);largest=Math.max(largest,style.radius);mean+=style.radius/types.length;
  });
  const first=frames[+slider.value];
  // omovi puts its point light half a box diagonal from the centre, which suits
  // Å-sized systems; in a dense reduced-unit box that is about 1σ in front of the
  // nearest face, so most spheres are lit from the side (half-moon shading). Use a
  // key light just above and left of the camera instead, without distance fall-off,
  // so boxes of every size are lit alike.
  const light=visualizer.pointLight;
  light.distance=0;light.decay=0;light.intensity=2.5;
  visualizer.updatePointLightPosition=()=>{
    const camera=visualizer.camera;camera.updateMatrixWorld();
    const right=new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld,0);
    const up=new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld,1);
    light.position.copy(camera.position).addScaledVector(up.sub(right),.3*camera.position.distanceTo(visualizer.getCameraTarget()));
  };
  visualizer.initPostProcessing({ssao:{enabled:true,radius:4*mean,intensity:5}});
  const particles=new lib.Particles(types.length);
  types.forEach((t,i)=>particles.add(...first.x[i].map((v,d)=>v-.5*first.box[d]),i,t));
  visualizer.add(particles);
  const pairs=data.bonds||(data.molecular?Array.from({length:Math.floor(types.length/2)},(_,k)=>[2*k,2*k+1]):[]);
  const bonds=pairs.length?new lib.Bonds(pairs.length):null;
  if(bonds)visualizer.add(bonds);
  // Bonds store a quarter of the given radius: pass 4x for a radius of 0.3 x the smallest atom.
  const bondRadius=4*.3*smallest;
  const box=new THREE.LineSegments(new THREE.BufferGeometry(),new THREE.LineBasicMaterial());
  // omovi draws opaque: paint the scene with the page colour, and box edges in a
  // grey that stands out from it. Follow light/dark theme switches.
  function matchPage(){
    const [r,g,b]=pageBackground(root.parentElement||root);
    visualizer.scene.background=new THREE.Color().setRGB(r/255,g/255,b/255,THREE.SRGBColorSpace);
    const light=.2126*r+.7152*g+.0722*b>128;
    box.material.color.setRGB(...(light?[.45,.5,.55]:[.62,.68,.73]),THREE.SRGBColorSpace);
    visualizer.forceRender=true;
  }
  matchPage();
  const theme=window.matchMedia?.('(prefers-color-scheme: dark)');
  theme?.addEventListener('change',matchPage);
  box.geometry.setAttribute('position',new THREE.Float32BufferAttribute(new Float32Array(72),3));
  box.frustumCulled=false;
  if(data.boundary!=='open')visualizer.scene.add(box);
  const aspect=()=>container.clientWidth/container.clientHeight||1.2;
  // A 35° field of view (omovi: 60°) keeps the perspective mild, so that a box
  // that fits in the view still fills it and spheres near the edges stay round.
  visualizer.perspectiveCamera.fov=35;
  // At zoom 1 all eight box corners are in view, with a 5% margin. Smaller zooms
  // step back; larger ones narrow the field of view, so the camera never moves
  // into the box.
  function resetCamera(){
    const lengths=frames[+slider.value].box,longest=lengths.indexOf(Math.max(...lengths));
    let direction=initial[2];
    if(projection.value==='orthographic'){
      // Straight along y is the pole of omovi's y-up orbit: step off it by a hair,
      // on the side that makes the front view's up direction point up on screen.
      const [,up,toward]=frontView();
      direction=toward.map((v,i)=>v-1e-6*Math.sign(toward[1])*up[i]);
    }else if(Math.max(...lengths)>1.5*Math.min(...lengths)){
      // Elongated box, such as a slab with vacuum: look across the long axis.
      direction=longest===2?[1,.35,0]:[0,.35,1];
      const norm=Math.hypot(...direction);direction=direction.map(v=>v/norm);
    }
    const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    const across=Math.hypot(direction[0],direction[2]),right=[direction[2]/across,0,-direction[0]/across];
    const up=[direction[1]*right[2]-direction[2]*right[1],direction[2]*right[0]-direction[0]*right[2],direction[0]*right[1]-direction[1]*right[0]];
    const camera=visualizer.perspectiveCamera,tan=Math.tan(Math.PI*camera.fov/360);
    let distance=0;
    for(let c=0;c<8;c++){
      const corner=lengths.map((v,d)=>(((c>>d)&1)-.5)*v),depth=dot(corner,direction);
      distance=Math.max(distance,depth+Math.abs(dot(corner,right))/(tan*aspect()),depth+Math.abs(dot(corner,up))/tan);
    }
    const magnification=data.zoom||1;
    camera.zoom=Math.max(magnification,1);camera.updateProjectionMatrix();
    visualizer.setCameraTarget(new THREE.Vector3(0,0,0));
    visualizer.setCameraPosition(new THREE.Vector3(...direction.map(v=>v*1.05*distance/Math.min(magnification,1))));
    fitOrthographic();
  }
  // omovi's orthographic frustum is half the camera distance high, far closer than the
  // perspective view: fit the box corners instead, with the same margin and zoom.
  function fitOrthographic(){
    if(visualizer.camera!==visualizer.orthographicCamera)return;
    const camera=visualizer.orthographicCamera;camera.updateMatrixWorld();
    const right=new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld,0);
    const up=new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld,1);
    const lengths=frames[+slider.value].box;let across=0,high=0;
    for(let c=0;c<8;c++){
      const corner=new THREE.Vector3(...lengths.map((v,d)=>(((c>>d)&1)-.5)*v));
      across=Math.max(across,Math.abs(corner.dot(right)));high=Math.max(high,Math.abs(corner.dot(up)));
    }
    const half=1.05*Math.max(high+largest,(across+largest)/aspect())/(data.zoom||1);  // atoms on the faces stick out
    Object.assign(camera,{top:half,bottom:-half,left:-half*aspect(),right:half*aspect(),zoom:1});
    camera.updateProjectionMatrix();
    visualizer.forceRender=true;
  }
  function show(f){
    const half=f.box.map(v=>.5*v);
    for(let i=0;i<types.length;i++)for(let d=0;d<3;d++)particles.positions[3*i+d]=f.x[i][d]-half[d];
    particles.markNeedsUpdate();
    if(bonds){
      bonds.count=0;
      for(const [i,j] of pairs){
        if(j>=types.length||!f.x[i].every((v,d)=>Math.abs(v-f.x[j][d])<half[d]))continue;
        bonds.add(...f.x[i].map((v,d)=>v-half[d]),...f.x[j].map((v,d)=>v-half[d]),bondRadius);
      }
      bonds.markNeedsUpdate();
    }
    const edge=box.geometry.attributes.position.array;let k=0;
    for(let a=0;a<8;a++)for(let d=0;d<3;d++){const b=a^(1<<d);if(a<b){
      for(const c of [a,b])for(let e=0;e<3;e++)edge[k++]=(((c>>e)&1)-.5)*f.box[e];
    }}
    box.geometry.attributes.position.needsUpdate=true;
    visualizer.forceRender=true;
    wake();
  }
  function setProjection(value){
    visualizer.setOrthographic(value==='orthographic');
    fitOrthographic();
    visualizer.setSSAO(value!=='orthographic');  // omovi's ambient occlusion needs the perspective camera
  }
  // Pause rendering while the view is scrolled out of sight. omovi only resizes while
  // rendering, so a resize, a new frame or pointer input renders a few frames anyway,
  // in case the visibility report is late or missing (a stale image is stretched).
  let visible=true;
  const visibility=typeof IntersectionObserver==='function'?new IntersectionObserver(entries=>{
    for(const entry of entries){visible=entry.isIntersecting;visualizer.idle=!visible;}
  }):null;
  if(visibility)visibility.observe(container);
  function wake(){
    if(!visualizer.idle)return;
    visualizer.idle=false;
    requestAnimationFrame(()=>requestAnimationFrame(()=>{visualizer.idle=!visible;}));
  }
  const resized=new ResizeObserver(wake);resized.observe(container);
  for(const name of ['pointerdown','pointermove','wheel'])container.addEventListener(name,wake,{passive:true});
  omovi={show,setProjection,resetCamera,dispose(){if(visibility)visibility.disconnect();resized.disconnect();theme?.removeEventListener('change',matchPage);visualizer.dispose();}};
  root.querySelector('[data-zoom-out]').hidden=root.querySelector('[data-zoom-in]').hidden=true;  // omovi zooms with the wheel
  showLegend(name=>{const c=data.units?lib.AtomTypes?.[name]?.color:null;return c?[c.r,c.g,c.b]:PALETTE[kindOf(name)];});
  resetCamera();setProjection(projection.value);draw();
  note.textContent+=' Rendered with OMOVI.';
}
if(data.renderer!=='basic'){
  const fallback=error=>console.warn('OMOVI unavailable; using the built-in renderer.',error);
  try{loadOmovi().then(startOmovi).catch(fallback);}catch(error){fallback(error);}
}
return {
  setFrame(frame){
    frames[0]=frame;
    requestDraw();
  },
  setStatus,
  dispose(){disposed=true;stop();if(omovi)omovi.dispose();else renderer.dispose();resizeObserver.disconnect();resolution.removeEventListener('change',watchResolution);if(pendingDraw)cancelAnimationFrame(pendingDraw);}
};
}

// Sphere impostors: one quad per atom, with ray/sphere intersection and surface
// depth per pixel. Centre-sorted discs cannot resolve dense, nearly coplanar atoms.
function createMDSphereRenderer(canvas){
  const gl=canvas.getContext('webgl2',{alpha:true,antialias:true,depth:true});
  if(!gl)throw new Error('The MD viewer requires WebGL 2. Enable graphics acceleration in your notebook browser/editor.');
  const programs=[],buffers=[],arrays=[];
  function program(vertex,fragment){
    function shader(type,source){
      const s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);
      if(!gl.getShaderParameter(s,gl.COMPILE_STATUS)){const error=gl.getShaderInfoLog(s);gl.deleteShader(s);throw new Error(error);}
      return s;
    }
    const vs=shader(gl.VERTEX_SHADER,vertex),fs=shader(gl.FRAGMENT_SHADER,fragment),p=gl.createProgram();
    gl.attachShader(p,vs);gl.attachShader(p,fs);gl.linkProgram(p);gl.deleteShader(vs);gl.deleteShader(fs);
    if(!gl.getProgramParameter(p,gl.LINK_STATUS))throw new Error(gl.getProgramInfoLog(p));
    programs.push(p);return p;
  }
  const common=`
precision highp float;
uniform vec2 u_origin, u_scale, u_viewport;
uniform bool u_perspective;
uniform float u_radius;
`;
  const spheres=program(`#version 300 es
${common}
layout(location=0) in vec3 a_center;
layout(location=1) in float a_kind;
layout(location=2) in float a_radius;
flat out vec3 v_center;
flat out float v_kind;
flat out float v_radius;
void main(){
  // Project the enclosing cube to conservatively bound the sphere, including
  // the off-axis perspective silhouette (a scaled centre-disc is insufficient).
  vec2 lo=vec2(1e10),hi=vec2(-1e10);
  for(int i=0;i<8;i++){
    vec3 corner=a_center+a_radius*vec3((i&1)==0?-1.:1.,(i&2)==0?-1.:1.,(i&4)==0?-1.:1.);
    vec2 p=corner.xy*(u_perspective?2.5/(2.5-corner.z):1.);
    lo=min(lo,p);hi=max(hi,p);
  }
  // One physical pixel of padding prevents raster rounding clipping the rim.
  lo-=1./(u_viewport*u_scale);hi+=1./(u_viewport*u_scale);
  vec2 corners[6]=vec2[6](vec2(0,0),vec2(1,0),vec2(0,1),vec2(0,1),vec2(1,0),vec2(1,1));
  vec2 p=mix(lo,hi,corners[gl_VertexID]);
  gl_Position=vec4(2.*(u_origin+u_scale*p)-1.,0.,1.);
  v_center=a_center;v_kind=a_kind;v_radius=a_radius;
}`,`#version 300 es
${common}
flat in vec3 v_center;
flat in float v_kind;
flat in float v_radius;
out vec4 color;
void main(){
  vec2 p=(gl_FragCoord.xy/u_viewport-u_origin)/u_scale;
  vec3 origin=u_perspective?vec3(0,0,2.5):vec3(p,2.5);
  vec3 direction=u_perspective?normalize(vec3(p,-2.5)):vec3(0,0,-1);
  vec3 offset=origin-v_center;
  float b=dot(offset,direction);
  float discriminant=b*b-dot(offset,offset)+v_radius*v_radius;
  if(discriminant<0.)discard;
  float t=-b-sqrt(discriminant);
  if(t<0.)discard;
  vec3 hit=origin+t*direction;
  vec3 normal=normalize(hit-v_center);
  gl_FragDepth=(2.5-hit.z)/5.;
  vec3 light=normalize(vec3(-.5,.65,1.));
  vec3 halfDirection=normalize(light-direction);
  float diffuse=max(dot(normal,light),0.);
  float specular=pow(max(dot(normal,halfDirection),0.),40.);
  float kind=mod(v_kind,7.);
  vec3 base=kind<.5?vec3(.16,.57,.72):kind<1.5?vec3(.84,.48,.19):kind<2.5?vec3(.25,.65,.4):kind<3.5?vec3(.6,.4,.75):kind<4.5?vec3(.8,.3,.35):kind<5.5?vec3(.5,.55,.6):vec3(.9,.91,.92);
  color=vec4(base*(.26+.74*diffuse)+vec3(.65)*specular,1.);
}`);
  const lines=program(`#version 300 es
${common}
layout(location=0) in vec3 a_position;
out vec3 v_position;
void main(){
  float w=u_perspective?(2.5-a_position.z)/2.5:1.;
  vec2 xy=(2.*u_origin-1.)*w+2.*u_scale*a_position.xy;
  gl_Position=vec4(xy,0.,w);
  v_position=a_position;
}`,`#version 300 es
precision highp float;
in vec3 v_position;
uniform vec3 u_color;
out vec4 color;
void main(){gl_FragDepth=(2.5-v_position.z)/5.;color=vec4(u_color,1.);}
`);
  function geometry(instanced){
    const vao=gl.createVertexArray(),buffer=gl.createBuffer();arrays.push(vao);buffers.push(buffer);
    gl.bindVertexArray(vao);gl.bindBuffer(gl.ARRAY_BUFFER,buffer);
    gl.enableVertexAttribArray(0);gl.vertexAttribPointer(0,3,gl.FLOAT,false,instanced?20:12,0);
    if(instanced){
      gl.vertexAttribDivisor(0,1);gl.enableVertexAttribArray(1);
      gl.vertexAttribPointer(1,1,gl.FLOAT,false,20,12);gl.vertexAttribDivisor(1,1);
      gl.enableVertexAttribArray(2);gl.vertexAttribPointer(2,1,gl.FLOAT,false,20,16);gl.vertexAttribDivisor(2,1);
    }
    return {vao,buffer};
  }
  const atomGeometry=geometry(true),lineGeometry=geometry(false);
  const uniforms=new Map(programs.map(p=>[p,Object.fromEntries(['origin','scale','viewport','perspective','radius','color'].map(n=>[n,gl.getUniformLocation(p,'u_'+n)]))]));
  function configure(p,scene){
    gl.useProgram(p);const u=uniforms.get(p);
    gl.uniform2f(u.origin,(360+scene.pan[0])/720,(300-scene.pan[1])/600);
    gl.uniform2f(u.scale,300*scene.zoom/720,300*scene.zoom/600);
    gl.uniform2f(u.viewport,canvas.width,canvas.height);
    gl.uniform1i(u.perspective,scene.perspective);gl.uniform1f(u.radius,scene.radius);
  }
  function drawLines(vertices,color,scene){
    if(!vertices.length)return;
    configure(lines,scene);gl.uniform3fv(uniforms.get(lines).color,color);
    gl.bindVertexArray(lineGeometry.vao);gl.bindBuffer(gl.ARRAY_BUFFER,lineGeometry.buffer);
    gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(vertices.flat()),gl.DYNAMIC_DRAW);
    gl.drawArrays(gl.LINES,0,vertices.length);
  }
  return {
    draw(scene){
      gl.viewport(0,0,canvas.width,canvas.height);gl.enable(gl.DEPTH_TEST);gl.depthFunc(gl.LESS);
      gl.clearColor(0,0,0,0);gl.clearDepth(1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
      drawLines(scene.edges,[.66,.72,.77],scene);drawLines(scene.bonds,[.49,.55,.59],scene);
      configure(spheres,scene);gl.bindVertexArray(atomGeometry.vao);gl.bindBuffer(gl.ARRAY_BUFFER,atomGeometry.buffer);
      const atoms=new Float32Array(scene.positions.length*5);
      scene.positions.forEach((p,i)=>{atoms.set(p,5*i);atoms[5*i+3]=scene.kinds?.[i]??(scene.molecular?i%2:0);atoms[5*i+4]=scene.radii?.[i]??scene.radius;});
      gl.bufferData(gl.ARRAY_BUFFER,atoms,gl.DYNAMIC_DRAW);
      gl.drawArraysInstanced(gl.TRIANGLES,0,6,scene.positions.length);
      gl.bindVertexArray(null);
    },
    dispose(){arrays.forEach(v=>gl.deleteVertexArray(v));buffers.forEach(b=>gl.deleteBuffer(b));programs.forEach(p=>gl.deleteProgram(p));}
  };
}
