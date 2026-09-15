// Browser regression harness. Load after viewer.js with #check and #results.
// Uses real GPU pixels; no screenshot comparison or external JS dependencies.
document.getElementById('check').onclick=()=>{
  const report=document.getElementById('results');
  try{
    const canvas=document.createElement('canvas');canvas.width=canvas.height=160;
    const renderer=createMDSphereRenderer(canvas),gl=canvas.getContext('webgl2');
    const positions=[[-.045,0,0],[.045,0,0]];
    const scene={positions,edges:[],bonds:[],radius:.08,zoom:3,pan:[0,0],perspective:false,molecular:false};
    function pixels(s){
      renderer.draw(s);const out=new Uint8Array(canvas.width*canvas.height*4);
      gl.readPixels(0,0,canvas.width,canvas.height,gl.RGBA,gl.UNSIGNED_BYTE,out);return out;
    }
    const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
    const unit=a=>{const l=Math.sqrt(dot(a,a));return a.map(v=>v/l);};
    function reference(x,y,s){
      const p=[((x+.5)/160-.5)/(300*s.zoom/720),((y+.5)/160-.5)/(300*s.zoom/600)];
      const o=s.perspective?[0,0,2.5]:[...p,2.5],d=s.perspective?unit([...p,-2.5]):[0,0,-1];
      let nearest=Infinity,normal;
      for(const c of s.positions){
        const oc=o.map((v,i)=>v-c[i]),b=dot(oc,d),disc=b*b-dot(oc,oc)+s.radius*s.radius;
        if(disc<0)continue;
        const t=-b-Math.sqrt(disc);
        if(t<nearest){nearest=t;normal=unit(o.map((v,i)=>v+t*d[i]-c[i]));}
      }
      if(!normal)return [0,0,0,0];
      const light=unit([-.5,.65,1]),half=unit(light.map((v,i)=>v-d[i]));
      const diffuse=Math.max(dot(normal,light),0),specular=Math.max(dot(normal,half),0)**40;
      return [...[.16,.57,.72].map(v=>Math.round(255*Math.min(1,v*(.26+.74*diffuse)+.65*specular))),255];
    }
    let comparisons=0,maxError=0,orderDifferences=0;
    for(const perspective of [false,true]){
      scene.perspective=perspective;
      for(const tilt of [-.001,0,.001]){
        scene.positions=positions.map(([x,y,z])=>[x*Math.cos(tilt),y,-x*Math.sin(tilt)]);
        const a=pixels(scene),b=pixels({...scene,positions:[...scene.positions].reverse()});
        for(let i=0;i<a.length;i++)if(a[i]!==b[i])orderDifferences++;
        // Both pixels lie in overlapping projected discs. The visible sphere
        // must be decided by surface depth, not one global centre ordering.
        for(const x of [74,85]){
          const expected=reference(x,80,scene),offset=(80*160+x)*4;
          for(let k=0;k<4;k++)maxError=Math.max(maxError,Math.abs(a[offset+k]-expected[k]));
          comparisons++;
        }
      }
    }
    renderer.dispose();
    if(orderDifferences!==0||maxError>3)throw new Error(`order differences=${orderDifferences}, colour error=${maxError}`);
    report.textContent=`PASS: ${comparisons} sphere-surface checks, both projections and near-coplanar angles.\nReversing draw order: ${orderDifferences} changed channels. Maximum ray-reference error: ${maxError}/255.`;
  }catch(error){report.textContent='FAIL: '+error.message;throw error;}
};
