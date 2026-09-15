/* Shared 3D renderer for live widgets and self-contained saved players. */
function createMDView(el,data){
// Own the sizing element: widget hosts may manage/reset styles on their `el`.
// A fixed CSS width scales with editor zoom; max-width only shrinks narrow panes.
const root=document.createElement('div'),frames=data.frames;
root.dataset.mdView='';
root.style.cssText='display:block;box-sizing:border-box;width:720px;max-width:100%;min-width:0;font:14px system-ui;color:#233744';
el.replaceChildren(root);
root.innerHTML="\n<canvas tabindex=\"0\" width=\"720\" height=\"600\" style=\"display:block;width:100%;height:auto;aspect-ratio:6/5;background:#f4f7fa;touch-action:none\" aria-label=\"Molecular trajectory, drag to rotate\"></canvas>\n<div style=\"display:flex;gap:12px;align-items:center;margin:10px 0\"><button type=\"button\">Play</button><button type=\"button\" data-reset>Reset view</button><input aria-label=\"Trajectory frame\" type=\"range\" min=\"0\" value=\"0\" style=\"flex:1\"><output></output></div><div style=\"display:flex;align-items:center;gap:12px;flex-wrap:wrap\"><label>Projection <select aria-label=\"Projection\"><option value=\"orthographic\">Orthographic</option><option value=\"perspective\">Perspective</option></select></label><label style=\"display:flex;align-items:center;gap:8px\">Zoom <button type=\"button\" aria-label=\"Zoom out\" data-zoom-out>\u2212</button><input type=\"range\" aria-label=\"Zoom\" min=\"25\" max=\"300\" step=\"1\" value=\"100\" style=\"width:130px\"><button type=\"button\" aria-label=\"Zoom in\" data-zoom-in>+</button><span data-zoom-value>100%</span></label></div><p data-note></p>\n";
root.querySelector('[data-note]').textContent=data.note;
let state=data.live?'Running':'',pendingDraw=0;
const canvas=root.querySelector('canvas'),renderer=createMDSphereRenderer(canvas),button=root.querySelector('button'),slider=root.querySelector('input'),label=root.querySelector('output');
const projection=root.querySelector('select');projection.value=data.projection;
const zoomSlider=root.querySelector('[aria-label="Zoom"]'),zoomLabel=root.querySelector('[data-zoom-value]');
let zoom=data.zoom;zoomSlider.value=Math.round(100*zoom);zoomLabel.textContent=Math.round(100*zoom)+'%';
if(data.live){button.hidden=true;slider.hidden=true;canvas.setAttribute('aria-label','Live molecular simulation, drag to rotate');}
slider.max=frames.length-1;const extent=Math.max(...frames.flatMap(f=>f.box));
function normalize(p,box){return p.map((v,d)=>(v-.5*box[d])/extent+.5);}
const pointers=new Map();
let pan=[0,0],timer=null;
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
function cameraPosition(p){
  const centered=p.map(v=>v-.5);
  return rotation.map(row=>row.reduce((sum,v,i)=>sum+v*centered[i],0));
}
// Keep projection, particles and bonds in the same logical coordinate system.
// Rebuild the backing bitmap at the actual display resolution after editor zoom
// or pane resizing; changing its attributes must not change the CSS aspect ratio.
function draw(){
const bounds=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1;
const width=Math.max(1,Math.round(bounds.width*dpr));
const height=Math.max(1,Math.round(bounds.height*dpr));
if(canvas.width!==width||canvas.height!==height){canvas.width=width;canvas.height=height;}
const f=frames[+slider.value],edges=[],bonds=[];
const corners=Array.from({length:8},(_,a)=>cameraPosition(normalize([a&1,(a>>1)&1,(a>>2)&1].map((v,d)=>v*f.box[d]),f.box)));
for(let a=0;a<8;a++)for(let d=0;d<3;d++){const b=a^(1<<d);if(a<b)edges.push(corners[a],corners[b]);}
const positions=f.x.map(p=>cameraPosition(normalize(p,f.box)));
if(data.molecular)for(let i=0;i+1<positions.length;i+=2){
  if(f.x[i].every((v,d)=>Math.abs(v-f.x[i+1][d])<.5*f.box[d]))bonds.push(positions[i],positions[i+1]);
}
renderer.draw({positions,edges,bonds,radius:.5/extent,zoom,pan,perspective:projection.value==='perspective',molecular:data.molecular});
updateLabel(f);}
function updateLabel(f){
  label.textContent=(state?state+' · ':'')+'Step '+f.step+' · t* '+f.time.toFixed(2)
    +(data.live?' · T* '+f.temperature.toFixed(3)+' · P* '+f.pressure.toPrecision(4):'');
}
function stop(){clearInterval(timer);timer=null;button.textContent='Play';}
button.onclick=()=>{if(timer){stop();return;}button.textContent='Pause';timer=setInterval(()=>{if(!root.isConnected){stop();return;}slider.value=(+slider.value+1)%frames.length;draw();},100);};
projection.onchange=draw;slider.oninput=draw;
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
function setZoom(value){zoom=Math.max(.25,Math.min(3,value));zoomSlider.value=Math.round(100*zoom);zoomLabel.textContent=Math.round(100*zoom)+'%';draw();}
zoomSlider.oninput=()=>setZoom(+zoomSlider.value/100);
root.querySelector('[data-zoom-out]').onclick=()=>setZoom(zoom/1.2);
root.querySelector('[data-zoom-in]').onclick=()=>setZoom(zoom*1.2);
canvas.addEventListener('wheel',e=>{e.preventDefault();setZoom(zoom*Math.exp(e.deltaY*.001));},{passive:false});
root.querySelector('[data-reset]').onclick=()=>{rotation=initial.map(row=>row.slice());pan=[0,0];setZoom(data.zoom);};
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
  dispose(){stop();renderer.dispose();resizeObserver.disconnect();resolution.removeEventListener('change',watchResolution);if(pendingDraw)cancelAnimationFrame(pendingDraw);}
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
flat out vec3 v_center;
flat out float v_kind;
void main(){
  // Project the enclosing cube to conservatively bound the sphere, including
  // the off-axis perspective silhouette (a scaled centre-disc is insufficient).
  vec2 lo=vec2(1e10),hi=vec2(-1e10);
  for(int i=0;i<8;i++){
    vec3 corner=a_center+u_radius*vec3((i&1)==0?-1.:1.,(i&2)==0?-1.:1.,(i&4)==0?-1.:1.);
    vec2 p=corner.xy*(u_perspective?2.5/(2.5-corner.z):1.);
    lo=min(lo,p);hi=max(hi,p);
  }
  // One physical pixel of padding prevents raster rounding clipping the rim.
  lo-=1./(u_viewport*u_scale);hi+=1./(u_viewport*u_scale);
  vec2 corners[6]=vec2[6](vec2(0,0),vec2(1,0),vec2(0,1),vec2(0,1),vec2(1,0),vec2(1,1));
  vec2 p=mix(lo,hi,corners[gl_VertexID]);
  gl_Position=vec4(2.*(u_origin+u_scale*p)-1.,0.,1.);
  v_center=a_center;v_kind=a_kind;
}`,`#version 300 es
${common}
flat in vec3 v_center;
flat in float v_kind;
out vec4 color;
void main(){
  vec2 p=(gl_FragCoord.xy/u_viewport-u_origin)/u_scale;
  vec3 origin=u_perspective?vec3(0,0,2.5):vec3(p,2.5);
  vec3 direction=u_perspective?normalize(vec3(p,-2.5)):vec3(0,0,-1);
  vec3 offset=origin-v_center;
  float b=dot(offset,direction);
  float discriminant=b*b-dot(offset,offset)+u_radius*u_radius;
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
  vec3 base=v_kind>.5?vec3(.84,.48,.19):vec3(.16,.57,.72);
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
    gl.enableVertexAttribArray(0);gl.vertexAttribPointer(0,3,gl.FLOAT,false,instanced?16:12,0);
    if(instanced){
      gl.vertexAttribDivisor(0,1);gl.enableVertexAttribArray(1);
      gl.vertexAttribPointer(1,1,gl.FLOAT,false,16,12);gl.vertexAttribDivisor(1,1);
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
      const atoms=new Float32Array(scene.positions.length*4);
      scene.positions.forEach((p,i)=>{atoms.set(p,4*i);atoms[4*i+3]=scene.molecular?i%2:0;});
      gl.bufferData(gl.ARRAY_BUFFER,atoms,gl.DYNAMIC_DRAW);
      gl.drawArraysInstanced(gl.TRIANGLES,0,6,scene.positions.length);
      gl.bindVertexArray(null);
    },
    dispose(){arrays.forEach(v=>gl.deleteVertexArray(v));buffers.forEach(b=>gl.deleteBuffer(b));programs.forEach(p=>gl.deleteProgram(p));}
  };
}
