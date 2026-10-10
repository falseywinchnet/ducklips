import {initialDesign,compile,validateDesign,MATERIALS,N,P,wav,clamp} from './engine.mjs';
const $=id=>document.getElementById(id), svg=$('tube');
let design=initialDesign();try{const saved=JSON.parse(localStorage.getItem('ducklips-flat-tube-v1'));if(saved)design=validateDesign(saved);}catch{}
let compiled=compile(design), selected=1,tool='anchors',history=[],context,node,running=false,state=null,drag=null,selectedBlock=-1,take=0,saveTimer,editTimer;
let captureContext=null;
let wave=new Float32Array(512),waveIndex=0,modeFrame=0;
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const X=x=>60+x/compiled.length*880,Y=y=>225-y*geometryScale(), width=()=>compiled.length;
function geometryScale(){return Math.min(10500,155/Math.max(.001,...compiled.upper,...compiled.lower));}
function checkpoint(){history.push(JSON.stringify(design));if(history.length>80)history.shift();}
function notice(text){$('notice').textContent=text;}
function refresh(send=true){
 compiled=compile(design);if(captureContext)captureContext.changes.push({milliseconds:performance.now()-captureContext.started,design:JSON.parse(JSON.stringify(design))});if(send&&node)node.port.postMessage({type:'configure',compiled});
 clearTimeout(saveTimer);saveTimer=setTimeout(()=>localStorage.setItem('ducklips-flat-tube-v1',JSON.stringify(design)),200);
 draw();panel();
}
function edited(){clearTimeout(editTimer);editTimer=setTimeout(()=>refresh(),35);draw();}
function panel(){
 const b=design.bands[selected];$('band-name').textContent=b.name.slice(5);
 $('band-description').textContent=selected===1?'The second band. Its two supported surfaces respond to pressure from both directions.':selected===0?'The chamber behind the band. Returning pressure stays connected to the constriction.':'The enclosing material carries motion while its shape determines the air passage.';
 $('material').value=b.material;
 for(const k of ['mass','elasticity','support']){$(k).value=b[k];$(k+'-value').textContent=b[k].toFixed(2)+' ×';}
 $('pressure').value=design.pressure;$('pressure-value').textContent=design.pressure+' Pa';$('wander').value=design.wander;$('wander-value').textContent=(design.wander*100).toFixed(1)+'%';$('monitor-value').textContent=Math.round(+$('monitor').value*100)+'%';
 $('frequencies').innerHTML=compiled.patches[selected].frequencies.map(f=>`<span>${f.toFixed(1)} Hz</span>`).join('');
 $('bands').innerHTML=design.bands.map((b,i)=>`<button data-band="${i}" class="${i===selected?'selected':''}" aria-pressed="${i===selected}">${esc(b.name)}<span>${MATERIALS[b.material].name}</span></button>`).join('');
 $('undo').disabled=!history.length;
}
function draw(){
 let content='<defs><pattern id="grid" width="40" height="40" patternUnits="userSpaceOnUse"><path d="M 40 0 L 0 0 0 40" fill="none" stroke="#26352a" stroke-width=".5"/></pattern></defs><rect width="1000" height="450" fill="url(#grid)"/><line x1="40" x2="965" y1="225" y2="225" stroke="#40513f" stroke-dasharray="3 7"/>';
 let start=0;
 const top=[],bottom=[],actualTop=[],actualBottom=[];
 for(let i=0;i<N;i++){
   const motion=state?.displacement,mode=$('mode-view').checked?Math.sin(modeFrame*.04)*.0006*compiled.patches[Math.floor(i/P)].phi[(i%P)*P]:0;
   top.push([X(compiled.x[i]),Y(compiled.upper[i])]);bottom.push([X(compiled.x[i]),Y(-compiled.lower[i])]);
   actualTop.push([X(compiled.x[i]),Y(compiled.upper[i]+(motion?.[i]||0)+mode)]);actualBottom.push([X(compiled.x[i]),Y(-compiled.lower[i]-(motion?.[N+i]||0)-mode)]);
 }
 const path=points=>points.map(p=>p.join(',')).join(' ');
 const outline=[...actualTop,...actualBottom.slice().reverse()];content+=`<polygon points="${path(outline)}" fill="#1c3030" stroke="none"/>`;
 for(let i=0;i<N;i++){
   const p=((state?.pressure[i]||0)+(state?.pressure[N+i]||0))/2, strength=Math.min(.72,Math.abs(p)/900+.07),color=p<0?'#c195bd':'#75b9ce';
   content+=`<rect class="cell" x="${X(compiled.x[i]-compiled.dx[i]/2)}" y="${actualTop[i][1]}" width="${compiled.dx[i]/compiled.length*880+1}" height="${Math.max(1,actualBottom[i][1]-actualTop[i][1])}" fill="${color}" opacity="${strength}"/>`;
 }
 for(let b=0;b<6;b++){
   const band=design.bands[b],left=X(start),right=X(start+band.length);content+=`<rect x="${left}" y="26" width="${right-left}" height="398" fill="${b===selected?'#a7bb66':'#000'}" opacity="${b===selected?'.045':'0'}"/><line x1="${left}" x2="${left}" y1="55" y2="400" stroke="#36432f" stroke-dasharray="2 6"/><text x="${(left+right)/2}" y="32" fill="#92a58a" font-size="10" text-anchor="middle">${String(b+1).padStart(2,'0')}</text>`;
   const points=actualTop.slice(b*P,b*P+P),lower=actualBottom.slice(b*P,b*P+P),color=MATERIALS[band.material].color;
   content+=`<polyline points="${path(points)}" fill="none" stroke="${color}" stroke-width="${band.radius<.001?2:7}" stroke-linecap="round" stroke-linejoin="round"/><polyline points="${path(lower)}" fill="none" stroke="${color}" stroke-width="${band.radius<.001?2:7}" stroke-linecap="round" stroke-linejoin="round"/>`;
   if(b<5){content+=`<line x1="${points.at(-1)[0]}" y1="${points.at(-1)[1]}" x2="${actualTop[(b+1)*P][0]}" y2="${actualTop[(b+1)*P][1]}" stroke="${color}" stroke-width="7"/><line x1="${lower.at(-1)[0]}" y1="${lower.at(-1)[1]}" x2="${actualBottom[(b+1)*P][0]}" y2="${actualBottom[(b+1)*P][1]}" stroke="${color}" stroke-width="7"/>`;}
   for(const side of ['upper','lower']){
     const pos=band[side+'Anchor']??.5,ax=left+(right-left)*(.125+.75*pos),target=Y((side==='upper'?1:-1)*(side==='upper'?compiled.upper:compiled.lower)[b*P+Math.round(pos*(P-1))]),ay=target+(side==='upper'?-34:34);
     content+=`<line x1="${ax}" y1="${target}" x2="${ax}" y2="${ay}" stroke="#a1af7d" stroke-width="1"/><path class="anchor" data-band="${b}" data-side="${side}" role="button" tabindex="0" aria-label="${esc(band.name.slice(5))} ${side} anchor" d="M${ax} ${ay-8}l8 8-8 8-8-8Z" fill="${b===selected?'#d8ed96':'#697e5b'}" stroke="#142015" stroke-width="2"/><text x="${ax}" y="${ay+(side==='upper'?-16:24)}" text-anchor="middle" font-size="9" fill="#8c9b81">${side==='upper'?'U':'L'}${b+1}</text>`;
   }
   start+=band.length;
 }
 if($('links').checked){for(let i=0;i<N-1;i++)for(let side=0;side<2;side++){
   const a=side?actualBottom[i]:actualTop[i],b=side?actualTop[i+1]:actualBottom[i+1];content+=`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="#79b4c8" opacity=".3" stroke-width=".8"/>`;
 }for(let i=0;i<N-P;i++){content+=`<path d="M${actualTop[i][0]} ${actualTop[i][1]-7} Q${(actualTop[i][0]+actualTop[i+P][0])/2} ${Math.min(actualTop[i][1],actualTop[i+P][1])-35} ${actualTop[i+P][0]} ${actualTop[i+P][1]-7}" fill="none" stroke="#d1b489" opacity=".27" stroke-width=".8"/>`;}}
 design.blocks.forEach((b,i)=>{const x=60+b.x*880,w=b.width*880;content+=`<g class="volume-block" data-block="${i}" role="button" tabindex="0" aria-label="Volume block ${i+1}, ${b.volume} millilitres"><rect x="${x-w/2}" y="204" width="${w}" height="42" rx="5" fill="#b7c77e" fill-opacity=".15" stroke="${selectedBlock===i?'#f1f8bf':'#9caf71'}" stroke-dasharray="4 3"/><text x="${x}" y="229" fill="#e0ebad" font-size="11" text-anchor="middle">+${b.volume} mL</text></g>`;});
 content+='<text x="16" y="214" fill="#aabd84" font-size="9">PRESSURE</text><path d="M18 230h25m-6-5 6 5-6 5" stroke="#d8ed96" fill="none"/><text x="948" y="214" fill="#aabd84" font-size="9">SAMPLE</text><path d="M950 230h25m-6-5 6 5-6 5" stroke="#d8ed96" fill="none"/><text x="60" y="422" fill="#76876e" font-size="10">0 mm</text><text x="940" y="422" text-anchor="end" fill="#76876e" font-size="10">'+(compiled.length*1000).toFixed(0)+' mm · side cutaway</text>';
 if($('mode-view').checked)content+='<text x="500" y="65" text-anchor="middle" fill="#d8ed96" font-size="11">Mode 1 shape preview · illustrative amplitude, no audio excitation</text>';
 svg.innerHTML=content;
 const gap=state?.gap??2*Math.min(...compiled.r.slice(4,8));$('opening').innerHTML=(gap*1000).toFixed(2)+' <small>mm</small>';
}
function coords(e){const p=svg.createSVGPoint();p.x=e.clientX;p.y=e.clientY;return p.matrixTransform(svg.getScreenCTM().inverse());}
svg.addEventListener('pointerdown',e=>{
 const anchor=e.target.closest('.anchor'),block=e.target.closest('.volume-block'),point=coords(e);checkpoint();
 if(anchor){selected=+anchor.dataset.band;drag={type:'anchor',band:selected,side:anchor.dataset.side,x:point.x,y:point.y,initial:design.bands[selected][anchor.dataset.side],scale:geometryScale()};selectedBlock=-1;}
 else if(block){selectedBlock=+block.dataset.block;drag={type:'block',index:selectedBlock};}
 else if(tool==='blocks'&&point.x>=60&&point.x<=940&&point.y>50&&point.y<400){design.blocks.push({x:clamp((point.x-60)/880,0,1),volume:+$('block-volume').value,width:.12});selectedBlock=design.blocks.length-1;drag={type:'block',index:selectedBlock};refresh();}
 else {history.pop();return;}
 svg.setPointerCapture(e.pointerId);panel();
});
svg.addEventListener('pointermove',e=>{
 if(!drag)return;const point=coords(e);
 if(drag.type==='anchor'){
  const b=design.bands[drag.band],sign=drag.side==='upper'?-1:1;
  b[drag.side]=clamp(drag.initial+sign*(point.y-drag.y)/drag.scale,-b.radius+.000025,.02);
  const start=design.bands.slice(0,drag.band).reduce((s,b)=>s+b.length,0);
  b[drag.side+'Anchor']=clamp(((point.x-60)/880*width()-start)/b.length,0,1);
 }else design.blocks[drag.index].x=clamp((point.x-60)/880,0,1);
 edited();
});
const endDrag=()=>{if(drag){drag=null;clearTimeout(editTimer);refresh();}};svg.addEventListener('pointerup',endDrag);svg.addEventListener('pointercancel',endDrag);
svg.addEventListener('keydown',e=>{
 const anchor=e.target.closest('.anchor');if(anchor&&['ArrowUp','ArrowDown','ArrowLeft','ArrowRight'].includes(e.key)){
 e.preventDefault();checkpoint();const b=design.bands[+anchor.dataset.band],side=anchor.dataset.side;
 if(e.key==='ArrowUp'||e.key==='ArrowDown')b[side]=clamp(b[side]+(e.key==='ArrowUp'?1:-1)*(side==='upper'?1:-1)*.0001,-b.radius+.000025,.02);
 else b[side+'Anchor']=clamp((b[side+'Anchor']??.5)+(e.key==='ArrowRight'?.03:-.03),0,1);refresh();svg.querySelector(`[data-band="${anchor.dataset.band}"][data-side="${side}"]`)?.focus();
 }
 const block=e.target.closest('.volume-block');if(block)selectedBlock=+block.dataset.block;
 if((e.key==='Delete'||e.key==='Backspace')&&selectedBlock>=0){e.preventDefault();checkpoint();design.blocks.splice(selectedBlock,1);selectedBlock=-1;refresh();}
});
$('bands').onclick=e=>{const b=e.target.closest('button');if(b){selected=+b.dataset.band;panel();draw();}};
for(const mode of ['anchors','blocks'])$(mode).onclick=()=>{tool=mode;for(const m of ['anchors','blocks']){$(m).classList.toggle('active',m===mode);$(m).setAttribute('aria-pressed',m===mode);} $('hint').textContent=mode==='blocks'?'Click inside to add volume. Drag blocks; Delete removes one.':'Drag a diamond inward to tighten its band. Arrow keys also work.';};
for(const k of ['mass','elasticity','support','pressure','wander']){
 $(k).addEventListener('pointerdown',checkpoint);$(k).addEventListener('keydown',e=>{if(e.key.startsWith('Arrow'))checkpoint();});
 $(k).oninput=()=>{if(['pressure','wander'].includes(k))design[k]=+$(k).value;else design.bands[selected][k]=+$(k).value;refresh();};
}
$('material').onchange=()=>{checkpoint();design.bands[selected].material=$('material').value;refresh();};
$('monitor').oninput=()=>{if(node)node.port.postMessage({type:'level',value:+$('monitor').value});panel();};
$('links').onchange=draw;$('mode-view').onchange=draw;
$('undo').onclick=()=>{if(history.length){design=JSON.parse(history.pop());state=null;refresh();}};
$('reset').onclick=()=>{checkpoint();design=initialDesign();state=null;refresh(false);if(node)node.port.postMessage({type:'reset',compiled});};
function download(data,type,name){const url=URL.createObjectURL(new Blob([data],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);}
$('save').onclick=()=>download(JSON.stringify(design,null,2),'application/json','ducklips-tube.json');
$('load').onchange=async e=>{try{const file=e.target.files[0];if(!file)return;if(file.size>100000)throw Error('Design file is too large.');const value=validateDesign(JSON.parse(await file.text()));checkpoint();design=value;state=null;refresh();notice('Tube loaded.');}catch(e){notice(e.message);}finally{$('load').value='';}};
async function audio(){
 if(context){await context.resume();return;}
 context=new AudioContext({sampleRate:48000});
 try{
  if(context.audioWorklet&&globalThis.isSecureContext){
   await context.audioWorklet.addModule('./worklet.mjs');
   node=new AudioWorkletNode(context,'flat-tube',{numberOfInputs:0,numberOfOutputs:1,outputChannelCount:[1]});
  }else{
   const {createLanAudio}=await import('./lan-audio.mjs');node=createLanAudio(context);
  }
  node.connect(context.destination);
 }catch(error){await context.close();context=null;node=null;throw error;}
 node.port.postMessage({type:'configure',compiled});node.port.postMessage({type:'level',value:+$('monitor').value});
 node.port.onmessage=({data:d})=>{if(d.type==='state'){state=d;wave=d.scope;waveIndex=d.scopeIndex;$('pressure-reading').textContent=Math.round(d.drive)+' Pa inlet';if(d.fault){running=false;node.port.postMessage({type:'run',value:false});$('start').textContent='Restart pressure';$('status').textContent='Model stopped';notice(d.fault);}}
 if(d.type==='recording'){addTake(d.samples,d.rate,'Live tube',captureContext);captureContext=null;$('record').disabled=false;$('record').textContent='Record 2 seconds';}
 };
 await context.resume();
}
$('start').onclick=async()=>{try{await audio();if(state?.fault){node.port.postMessage({type:'reset',compiled});state=null;}running=!running;node.port.postMessage({type:'run',value:running});$('start').textContent=running?'Stop pressure':'Start pressure';$('status').textContent=running?'Running · outlet connected':'Pressure off · releasing';}catch(e){notice('Audio could not start: '+e.message);}};
$('record').onclick=async()=>{try{await audio();if(!running){running=true;node.port.postMessage({type:'run',value:true});$('start').textContent='Stop pressure';$('status').textContent='Running · recording outlet';}$('record').disabled=true;$('record').textContent='Recording…';captureContext={started:performance.now(),date:new Date().toISOString(),design:JSON.parse(JSON.stringify(design)),initialObservation:state,changes:[]};node.port.postMessage({type:'capture',seconds:2});}catch(e){notice(e.message);$('record').disabled=false;}};
function addTake(samples,rate,title,provenance){
 take++;const bytes=wav(samples,rate),blob=new Blob([bytes],{type:'audio/wav'}),url=URL.createObjectURL(blob);let peak=0;for(const x of samples)peak=Math.max(peak,Math.abs(x));
 const row=document.createElement('div');row.className='take';const label=document.createElement('span');label.textContent=`${String(take).padStart(2,'0')} · ${title}`;const player=document.createElement('audio');player.controls=true;player.src=url;const link=document.createElement('a');link.href=url;link.download=`ducklips-tube-${take}.wav`;link.textContent='WAV';const enhanced=document.createElement('button');enhanced.textContent='Depict enhanced';const receipt=document.createElement('small');receipt.textContent=`${(samples.length/rate).toFixed(2)} s · ${(rate/1000).toFixed(1)} kHz · peak ${peak.toFixed(4)}`;
 row.append(label,player,link,enhanced,receipt);if(take===1)$('recordings').textContent='';$('recordings').prepend(row);
 // Full wave, no decimation-by-selection: min/max bins preserve all extrema.
 wave=samples;waveIndex=0;
 enhanced.onclick=async()=>{enhanced.disabled=true;enhanced.textContent='Computing full field…';try{const response=await fetch('/api/enhanced',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({audio:await new Promise((resolve,reject)=>{const reader=new FileReader();reader.onload=()=>resolve(reader.result.split(',')[1]);reader.onerror=reject;reader.readAsDataURL(blob);}),provenance})});const result=await response.json();if(!response.ok)throw Error(result.error);const a=document.createElement('a');a.href=result.url;a.target='_blank';a.textContent='Full enhanced field ↗';enhanced.replaceWith(a);}catch(e){enhanced.disabled=false;enhanced.textContent='Retry enhanced';notice(e.message);}};
 // Always produce the paired full depiction for a recorded experiment.
 enhanced.click();
}
function scope(){const canvas=$('scope'),dpr=devicePixelRatio||1,w=canvas.clientWidth,h=112;if(canvas.width!==Math.round(w*dpr)){canvas.width=w*dpr;canvas.height=h*dpr;}const g=canvas.getContext('2d');g.setTransform(dpr,0,0,dpr,0,0);g.clearRect(0,0,w,h);g.strokeStyle='#34402e';g.beginPath();g.moveTo(0,h/2);g.lineTo(w,h/2);g.stroke();let peak=.002;for(const v of wave)peak=Math.max(peak,Math.abs(v));g.strokeStyle='#bfd693';g.lineWidth=1;g.beginPath();for(let x=0;x<w;x++){let lo=0,hi=0;const a=Math.floor(x*wave.length/w),b=Math.max(a+1,Math.floor((x+1)*wave.length/w));for(let j=a;j<b;j++){const v=wave[(j+waveIndex)%wave.length];lo=Math.min(lo,v);hi=Math.max(hi,v);}g.moveTo(x,h/2-hi/peak*45);g.lineTo(x,h/2-lo/peak*45);}g.stroke();g.fillStyle='#869679';g.font='10px sans-serif';g.fillText(`±${peak.toFixed(4)} · display autoscale`,12,16);}
let last=0;function animate(time){if(time-last>40){last=time;modeFrame++;if(running||(state?.drive>0.1)||(state?.velocity.some(v=>Math.abs(v)>1e-6))||$('mode-view').checked)draw();scope();}requestAnimationFrame(animate);}refresh(false);requestAnimationFrame(animate);
// Benchmark the actual core away from the UI thread. No device or spectral work included.
const worker=new Worker('./render-worker.mjs',{type:'module'});worker.onmessage=({data:d})=>{if(d.error){notice(d.error);return;}$('speed').textContent=d.realtime.toFixed(1)+'× core realtime';if(d.realtime<10)notice('Core is below the 10× target on this browser. Geometry compilation and drawing are measured separately.');worker.terminate();};worker.postMessage({compiled:compile({...design,pressure:0,wander:0}),seconds:1,rate:24000});
