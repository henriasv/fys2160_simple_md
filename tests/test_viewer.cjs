// Canvas event tests; run with: node md-simulator/tests/test_viewer.cjs
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
let images = [], disconnected = false;
const context = {
  clearRect(){ images=[]; }, setTransform(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){},
  createRadialGradient(){return {addColorStop(){}};}, arc(){}, fill(){},
  drawImage(sprite,x,y,w,h){images.push([x+w/2,y+h/2,w]);}
};
const canvas = {
  width:720,height:600,setAttribute(){},setPointerCapture(){},
  getContext(){return context;},getBoundingClientRect(){return {width:720,height:600};},
  addEventListener(name,fn){this['on'+name]=fn;}
};
const controls = new Map([['canvas',canvas]]);
const root = {style:{},dataset:{},isConnected:true,
  querySelector(selector){if(!controls.has(selector))controls.set(selector,{value:0});return controls.get(selector);}
};
const host = {style:{cssText:'display:contents'},replaceChildren(child){this.child=child;}};
const sandbox = {
  document:{createElement(tag){return tag==='div'?root:{getContext(){return context;}};}},
  window:{devicePixelRatio:2,matchMedia(){return {addEventListener(){},removeEventListener(){}};}},
  ResizeObserver:class {observe(){} disconnect(){disconnected=true;}},
  requestAnimationFrame(){return 1;},cancelAnimationFrame(){},setInterval,clearInterval
};
vm.createContext(sandbox);
vm.runInContext(fs.readFileSync(__dirname+'/../src/fys2160_md/viewer.js','utf8'),sandbox);
const frame={box:[10,10,10],x:[[5,5,5],[7,6,5]],step:0,time:0,temperature:2,pressure:1};
const view=sandbox.createMDView(host,{frames:[frame],projection:'orthographic',zoom:1,live:true,note:''});
assert.equal(host.child,root);
assert.equal(host.style.cssText,'display:contents'); // Renderer leaves host layout alone.
assert.equal(canvas.width,1440); // Retina backing bitmap.
const original=images.map(p=>p.slice());
const event=(id,x,y,button=0,shiftKey=false)=>({pointerId:id,clientX:x,clientY:y,button,shiftKey,preventDefault(){}});
function drag(button,shift=false){
  canvas.onpointerdown(event(1,100,100,button,shift));
  canvas.onpointermove(event(1,140,120,button,shift));
  canvas.onpointerup(event(1,140,120,button,shift));
}
function reset(){controls.get('[data-reset]').onclick();assert.deepEqual(images,original);}
for(const [button,shift] of [[2,false],[0,true]]){
  drag(button,shift);
  assert.equal(images[0][0]-original[0][0],40);
  assert.equal(images[0][1]-original[0][1],20);
  assert.equal(images[0][2],original[0][2]);
  reset();
}
drag(0);
assert.notDeepEqual(images,original); // Primary drag rotates off-centre atoms.
reset();
canvas.onpointerdown(event(1,100,100));canvas.onpointerdown(event(2,200,100));
canvas.onpointermove(event(1,140,120));canvas.onpointermove(event(2,240,120));
assert.equal(images[0][0]-original[0][0],40);
assert.equal(images[0][1]-original[0][1],20);
canvas.onpointercancel(event(1,140,120));canvas.onpointerup(event(2,240,120));
reset();
canvas.onwheel({deltaY:100,preventDefault(){}});
assert.ok(images[0][2]>original[0][2]); // Requested inversion: positive scroll zooms in.
reset();
view.dispose();assert.ok(disconnected);
const expanded={...frame,box:[20,20,20],x:frame.x.map(p=>p.map(v=>2*v))};
const second=sandbox.createMDView(host,{frames:[expanded],projection:'orthographic',zoom:1,live:true,note:''});
assert.equal(images[0][2],original[0][2]/2); // Same N, 8x lower density: half the projected radius.
assert.equal(original[0][2]*60/64,30); // Diameter sigma=1 at 300/10 pixels per sigma.
second.dispose();
console.log('Viewer controls passed: rotation, secondary/Shift/two-touch pan, inverted zoom, reset, sizing and cleanup.');
