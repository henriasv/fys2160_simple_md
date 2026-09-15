/* Shared 3D renderer for live widgets and self-contained saved players. */
function createMDView(el,data){
const root=el,frames=data.frames;
root.style.cssText='width:100%;max-width:52em;min-width:0;font-family:system-ui;color:#233744';
root.innerHTML="\n<canvas tabindex=\"0\" width=\"720\" height=\"600\" style=\"display:block;width:100%;height:auto;aspect-ratio:6/5;background:#f4f7fa;touch-action:none\" aria-label=\"Molecular trajectory, drag to rotate\"></canvas>\n<div style=\"display:flex;gap:12px;align-items:center;margin:10px 0\"><button type=\"button\">Play</button><button type=\"button\" data-reset>Reset view</button><input aria-label=\"Trajectory frame\" type=\"range\" min=\"0\" value=\"0\" style=\"flex:1\"><output></output></div><div style=\"display:flex;align-items:center;gap:12px;flex-wrap:wrap\"><label>Projection <select aria-label=\"Projection\"><option value=\"orthographic\">Orthographic</option><option value=\"perspective\">Perspective</option></select></label><label style=\"display:flex;align-items:center;gap:8px\">Zoom <button type=\"button\" aria-label=\"Zoom out\" data-zoom-out>\u2212</button><input type=\"range\" aria-label=\"Zoom\" min=\"25\" max=\"300\" step=\"1\" value=\"100\" style=\"width:130px\"><button type=\"button\" aria-label=\"Zoom in\" data-zoom-in>+</button><span data-zoom-value>100%</span></label></div><p data-note></p>\n";
root.querySelector('[data-note]').textContent=data.note;
let state=data.live?'Running':'',pendingDraw=0;
const canvas=root.querySelector('canvas'),ctx=canvas.getContext('2d'),button=root.querySelector('button'),slider=root.querySelector('input'),label=root.querySelector('output');
const projection=root.querySelector('select');projection.value=data.projection;
const zoomSlider=root.querySelector('[aria-label="Zoom"]'),zoomLabel=root.querySelector('[data-zoom-value]');
let zoom=data.zoom;zoomSlider.value=Math.round(100*zoom);zoomLabel.textContent=Math.round(100*zoom)+'%';
if(data.live){button.hidden=true;slider.hidden=true;canvas.setAttribute('aria-label','Live molecular simulation, drag to rotate');}
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
// Keep projection, particles and bonds in the same logical coordinate system.
// Rebuild the backing bitmap at the actual display resolution after editor zoom
// or pane resizing; changing its attributes must not change the CSS aspect ratio.
function draw(){
const bounds=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1;
const width=Math.max(1,Math.round(bounds.width*dpr));
const height=Math.max(1,Math.round(bounds.height*dpr));
if(canvas.width!==width||canvas.height!==height){canvas.width=width;canvas.height=height;}
ctx.setTransform(width/720,0,0,height/600,0,0);
const f=frames[+slider.value];ctx.clearRect(0,0,720,600);
for(let a=0;a<8;a++)for(let d=0;d<3;d++){const b=a^(1<<d);if(a<b)line(project(normalize([a&1,(a>>1)&1,(a>>2)&1].map((v,d)=>v*f.box[d]),f.box)),project(normalize([b&1,(b>>1)&1,(b>>2)&1].map((v,d)=>v*f.box[d]),f.box)),'#a8b8c4',1);}
const positions=f.x.map(p=>p.map((v,d)=>v/f.box[d]));const points=f.x.map(p=>project(normalize(p,f.box)));
if(data.molecular)for(let i=0;i+1<points.length;i+=2){if(positions[i].every((v,d)=>Math.abs(v-positions[i+1][d])<.5))line(points[i],points[i+1],'#7d8b96',2);}
points.map((p,i)=>({p,i})).sort((a,b)=>a.p[2]-b.p[2]).forEach(({p,i})=>{const radius=4.5*zoom*p[3];ctx.drawImage(sprites[data.molecular&&i%2?1:0],p[0]-radius,p[1]-radius,2*radius,2*radius);});
updateLabel(f);}
function updateLabel(f){
  label.textContent=(state?state+' · ':'')+'Step '+f.step+' · t* '+f.time.toFixed(2)
    +(data.live?' · T* '+f.temperature.toFixed(3)+' · P* '+f.pressure.toPrecision(4):'');
}
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
return {
  setFrame(frame){
    frames[0]=frame;
    requestDraw();
  },
  setStatus(status){state=status.charAt(0).toUpperCase()+status.slice(1);updateLabel(frames[0]);},
  dispose(){stop();resizeObserver.disconnect();resolution.removeEventListener('change',watchResolution);if(pendingDraw)cancelAnimationFrame(pendingDraw);}
};
}
