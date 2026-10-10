import {systemDesign} from './system-presets.mjs';
import {compileSystem,validateSystem} from './system.mjs';
import {referenceDesign} from './presets.mjs';
import {initialDesign,compile,validateDesign,MATERIALS,C,modeShape,wav,clamp} from './engine.mjs';
const $=id=>document.getElementById(id), svg=$('tube');
let system=systemDesign('open');try{const saved=JSON.parse(localStorage.getItem('ducklips-connected-v4'));if(saved)system=validateSystem(saved);}catch{}
let stage=0,design=system.tubes?.[stage].design??system,whole=compileSystem(system),compiled=whole.tubes?.[stage].compiled??whole,N=compiled.n,selected=compiled.cordIndex,history=[],context,node,running=false,state=null,systemState=null,drag=null,take=0,saveTimer,editTimer;
$('monitor').value=system.version===4?8:.25;
let captureContext=null,modeIndex=0,previewShape=modeShape(compiled.sheet,0);
let wave=new Float32Array(512),waveIndex=0,modeFrame=0;
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const X=x=>60+x/compiled.length*880,Y=y=>225-y*geometryScale(), width=()=>compiled.length;
function geometryScale(){return Math.min(10500,155/Math.max(.001,...compiled.upper,...compiled.lower));}
function checkpoint(){history.push(JSON.stringify(system));if(history.length>80)history.shift();}
function notice(text){$('notice').textContent=text;}
function refresh(send=true){
 whole=compileSystem(system);stage=Math.min(stage,(system.tubes?.length??1)-1);design=system.tubes?.[stage].design??system;compiled=whole.tubes?.[stage].compiled??whole;state=systemState?.tubes?.[stage]??systemState;N=compiled.n;selected=Math.min(selected,design.bands.length-1);modeIndex=Math.min(modeIndex,compiled.sheet.modes.length-1);previewShape=modeShape(compiled.sheet,modeIndex);if(captureContext)captureContext.changes.push({milliseconds:performance.now()-captureContext.started,design:JSON.parse(JSON.stringify(system))});if(send&&node)node.port.postMessage({type:'configure',compiled:whole});
 clearTimeout(saveTimer);saveTimer=setTimeout(()=>localStorage.setItem('ducklips-connected-v4',JSON.stringify(system)),200);
 draw();panel();network();
}
function edited(){clearTimeout(editTimer);editTimer=setTimeout(()=>refresh(),35);draw();}
function panel(){
 const b=design.bands[selected];$('band-name').textContent=b.name.slice(5);
 $('band-description').textContent=(whole.kind!=='network'||stage===whole.source)&&selected===compiled.cordIndex?'The constricting elastic cord. Feedback arrives through connected tissue and air on both sides.':design.bands[selected].role==='adapter'?'A light elastic adapter with little support, open to the main tube diameter.':design.bands[selected].role==='soft'?'Heavy soft tissue with little external support; it responds to returning pressure.':'This section shapes returning pressure and carries material motion.';
 $('material').value=b.material;$('preset').value=system.version===4?system.arrangement:system.version===3?'below-mouth':'reference';$('bands').style.setProperty('--sections',design.bands.length);
 $('body-elasticity').value=system.elasticityScale??1;$('body-elasticity-value').textContent=(system.elasticityScale??1).toFixed(2)+' ×';
 $('elasticity-variation').value=system.elasticityVariation??0;$('elasticity-variation-value').textContent='±'+Math.round((system.elasticityVariation??0)*100)+'% maximum';
 for(const k of ['mass','elasticity','support']){$(k).value=b[k];$(k+'-value').textContent=b[k].toFixed(2)+' ×';}
 $('pressure').value=system.pressure;$('pressure-value').textContent=system.pressure+' Pa';$('wander').value=system.wander;$('wander-value').textContent=(system.wander*100).toFixed(1)+'%';$('monitor-value').textContent=(+$('monitor').value).toFixed(2)+' ×';
 $('mode-select').innerHTML=compiled.sheet.modes.map((m,i)=>`<option value="${i}">${String(i+1).padStart(3,'0')} · ${m.frequency.toFixed(1)} Hz · ${m.k===0?'ring motion':m.k+' around circumference'}</option>`).join('');$('mode-select').value=modeIndex;
 $('mode-summary').textContent=`${compiled.sheet.modes.length} nominal modes · ${compiled.sheet.modes[0].frequency.toFixed(1)}–${compiled.sheet.modes.at(-1).frequency.toFixed(1)} Hz`;
 $('bands').innerHTML=design.bands.map((b,i)=>`<button data-band="${i}" class="${i===selected?'selected':''}" aria-pressed="${i===selected}">${esc(b.name)}<span>${MATERIALS[b.material].name}</span></button>`).join('');
 $('undo').disabled=!history.length;
}
function network(){
 const connected=system.version===4;$('network').hidden=!connected;if(!connected)return;
 const buttons=system.tubes.map((t,i)=>`<button data-stage="${i}" class="${i===stage?'selected':''}" aria-pressed="${i===stage}">${esc(t.name)}<small>${t.design.bands.length} sections · ${whole.tubes[i].compiled.sheet.modes.length} modes</small></button>`);
 $('stages').innerHTML=buttons.join('');
 $('stage-title').textContent=system.tubes[stage].name;
 $('network-summary').textContent=`${whole.tubes.reduce((n,t)=>n+t.compiled.sheet.modes.length,0)} sheet modes · pressure and tissue feedback travel both ways · sample at the lips`;
 $('route').innerHTML=system.connections.map(e=>`<span>${esc(system.tubes.find(t=>t.id===e.a.tube).name)} ${e.a.section!==undefined?'section '+(e.a.section+1):''} ↔ ${esc(system.tubes.find(t=>t.id===e.b.tube).name)}</span>`).join('');
}
$('stages').onclick=e=>{const button=e.target.closest('[data-stage]');if(!button)return;stage=+button.dataset.stage;selected=0;modeIndex=0;design=system.tubes[stage].design;compiled=whole.tubes[stage].compiled;N=compiled.n;state=systemState?.tubes?.[stage]??null;previewShape=modeShape(compiled.sheet,0);draw();panel();network();};
function draw(){
 let content='<defs><pattern id="grid" width="40" height="40" patternUnits="userSpaceOnUse"><path d="M 40 0 L 0 0 0 40" fill="none" stroke="#26352a" stroke-width=".5"/></pattern></defs><rect width="1000" height="450" fill="url(#grid)"/><line x1="40" x2="965" y1="225" y2="225" stroke="#40513f" stroke-dasharray="3 7"/>';
 let start=0;
 const top=[],bottom=[],actualTop=[],actualBottom=[];
 for(let i=0;i<N;i++){
   const motion=state?.displacement,amplitude=$('mode-view').checked?Math.sin(modeFrame*.04)*.0006:0,mode=amplitude*previewShape[i],lowerMode=amplitude*previewShape[C/2*N+i];
   top.push([X(compiled.x[i]),Y(compiled.upper[i])]);bottom.push([X(compiled.x[i]),Y(-compiled.lower[i])]);
   actualTop.push([X(compiled.x[i]),Y(compiled.upper[i]+(motion?.[i]||0)+mode)]);actualBottom.push([X(compiled.x[i]),Y(-compiled.lower[i]-(motion?.[N+i]||0)-lowerMode)]);
 }
 const path=points=>points.map(p=>p.join(',')).join(' ');
 const outline=[...actualTop,...actualBottom.slice().reverse()];content+=`<polygon points="${path(outline)}" fill="#1c3030" stroke="none"/>`;
 for(let i=0;i<N;i++){
   const p=((state?.pressure[i]||0)+(state?.pressure[N+i]||0))/2, strength=Math.min(.72,Math.abs(p)/900+.07),color=p<0?'#c195bd':'#75b9ce';
   content+=`<rect class="cell" x="${X(compiled.x[i]-compiled.dx[i]/2)}" y="${actualTop[i][1]}" width="${compiled.dx[i]/compiled.length*880+1}" height="${Math.max(1,actualBottom[i][1]-actualTop[i][1])}" fill="${color}" opacity="${strength}"/>`;
 }
 for(let b=0;b<design.bands.length;b++){
   const band=design.bands[b],first=compiled.starts[b],count=compiled.counts[b],left=X(start),right=X(start+band.length);content+=`<rect x="${left}" y="26" width="${right-left}" height="398" fill="${b===selected?'#a7bb66':'#000'}" opacity="${b===selected?'.045':'0'}"/><line x1="${left}" x2="${left}" y1="55" y2="400" stroke="#36432f" stroke-dasharray="2 6"/><text x="${(left+right)/2}" y="32" fill="#92a58a" font-size="10" text-anchor="middle">${String(b+1).padStart(2,'0')}</text>`;
   const points=actualTop.slice(first,first+count),lower=actualBottom.slice(first,first+count),color=MATERIALS[band.material].color;
   content+=`<polyline points="${path(points)}" fill="none" stroke="${color}" stroke-width="${band.radius<.001?2:7}" stroke-linecap="round" stroke-linejoin="round"/><polyline points="${path(lower)}" fill="none" stroke="${color}" stroke-width="${band.radius<.001?2:7}" stroke-linecap="round" stroke-linejoin="round"/>`;
   if(b<design.bands.length-1){content+=`<line x1="${points.at(-1)[0]}" y1="${points.at(-1)[1]}" x2="${actualTop[first+count][0]}" y2="${actualTop[first+count][1]}" stroke="${color}" stroke-width="7"/><line x1="${lower.at(-1)[0]}" y1="${lower.at(-1)[1]}" x2="${actualBottom[first+count][0]}" y2="${actualBottom[first+count][1]}" stroke="${color}" stroke-width="7"/>`;}
   const pos=band.anchor,ax=left+(right-left)*(.5/count+(1-1/count)*pos),target=Y(compiled.r[first+Math.round(pos*(count-1))]),ay=target-34;
   content+=`<line x1="${ax}" y1="${target}" x2="${ax}" y2="${ay}" stroke="#a1af7d" stroke-width="1"/><path class="anchor" data-band="${b}" role="button" tabindex="0" aria-label="${esc(band.name.slice(5))} ring support" d="M${ax} ${ay-8}l8 8-8 8-8-8Z" fill="${b===selected?'#d8ed96':'#697e5b'}" stroke="#142015" stroke-width="2"/><text x="${ax}" y="${ay-16}" text-anchor="middle" font-size="9" fill="#8c9b81">${b+1}</text>`;
   start+=band.length;
 }
 if($('links').checked){for(let i=0;i<N-1;i++)for(let side=0;side<2;side++){
   const a=side?actualBottom[i]:actualTop[i],b=side?actualTop[i+1]:actualBottom[i+1];content+=`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="#79b4c8" opacity=".3" stroke-width=".8"/>`;
 }for(let i=0;i<N-1;i++){content+=`<path d="M${actualTop[i][0]} ${actualTop[i][1]-7} Q${(actualTop[i][0]+actualTop[i+1][0])/2} ${Math.min(actualTop[i][1],actualTop[i+1][1])-35} ${actualTop[i+1][0]} ${actualTop[i+1][1]-7}" fill="none" stroke="#d1b489" opacity=".27" stroke-width=".8"/>`;}}
 content+='<text x="16" y="214" fill="#aabd84" font-size="9">'+(whole.kind==='network'&&stage!==whole.source?'RETURN':'PRESSURE')+'</text><path d="M18 230h25m-6-5 6 5-6 5" stroke="#d8ed96" fill="none"/><text x="948" y="214" fill="#aabd84" font-size="9">'+(whole.kind==='network'&&stage!==whole.output?(system.tubes[stage].end==='closed'?'CLOSED':'NEXT'):'SAMPLE')+'</text><path d="M950 230h25m-6-5 6 5-6 5" stroke="#d8ed96" fill="none"/><text x="60" y="422" fill="#76876e" font-size="10">0 mm</text><text x="940" y="422" text-anchor="end" fill="#76876e" font-size="10">'+(compiled.length*1000).toFixed(0)+' mm · side cutaway</text>';
 if(whole.kind==='network'){
  for(const link of whole.connections)for(const [port,other] of [[link.a,link.b],[link.b,link.a]])if(port.tube===stage&&!port.end){
   const x=X(compiled.x[port.index]),y=actualBottom[port.index][1];content+=`<path d="M${x} ${y}v32" stroke="#c49ad7" stroke-width="2"/><circle cx="${x}" cy="${y}" r="4" fill="#c49ad7"/><text x="${x}" y="${y+46}" fill="#c49ad7" font-size="10" text-anchor="middle">↔ ${esc(system.tubes[other.tube].name)}</text>`;
  }
  if(system.tubes[stage].end==='closed'){const a=actualTop.at(-1),b=actualBottom.at(-1);content+=`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="#8da590" stroke-width="5"/>`;}
 }
 if($('mode-view').checked)content+=`<text x="500" y="65" text-anchor="middle" fill="#d8ed96" font-size="11">Tube mode ${modeIndex+1} · ${compiled.sheet.modes[modeIndex].frequency.toFixed(1)} Hz · shape preview</text>`;
 svg.innerHTML=content;drawSheet();
 const body=whole.tubes?.[whole.source].compiled??whole;const gap=systemState?.gap??2*Math.min(...body.r.slice(body.starts[body.cordIndex],body.starts[body.cordIndex]+body.counts[body.cordIndex]));$('opening').innerHTML=(gap*1000).toFixed(2)+' <small>mm</small>';
}
function drawSheet(){
 const shown=$('mode-view').checked;$('sheet-view').hidden=!shown;if(!shown)return;
 const canvas=$('sheet-map'),g=canvas.getContext('2d');g.clearRect(0,0,960,180);
 for(let t=0;t<C;t++)for(let i=0;i<N;i++){
  const value=previewShape[t*N+i],strength=Math.abs(value),x=60+(compiled.x[i]-compiled.dx[i]/2)/compiled.length*880;
  g.fillStyle=value>=0?`rgb(${Math.round(35+190*strength)},${Math.round(50+125*strength)},${Math.round(45+65*strength)})`:`rgb(${Math.round(25+60*strength)},${Math.round(45+110*strength)},${Math.round(50+165*strength)})`;
  g.fillRect(x,20+t*17.5,compiled.dx[i]/compiled.length*880+.3,17.5);
 }
 g.strokeStyle='#d8ed96';g.setLineDash([4,4]);g.beginPath();g.moveTo(60,18);g.lineTo(940,18);g.moveTo(60,162);g.lineTo(940,162);g.stroke();g.setLineDash([]);
 g.fillStyle='#aab99f';g.font='11px sans-serif';g.fillText('0°',15,28);g.fillText('360°',6,163);
 $('sheet-caption').textContent=`Mode ${modeIndex+1} of ${compiled.sheet.modes.length} · ${compiled.sheet.modes[modeIndex].frequency.toFixed(1)} Hz. The dashed edges join; colors show the signed radial deformation over the complete sheet.`;
}
function coords(e){const p=svg.createSVGPoint();p.x=e.clientX;p.y=e.clientY;return p.matrixTransform(svg.getScreenCTM().inverse());}
svg.addEventListener('pointerdown',e=>{
 const anchor=e.target.closest('.anchor'),point=coords(e);checkpoint();
 if(anchor){selected=+anchor.dataset.band;drag={type:'anchor',band:selected,x:point.x,y:point.y,initial:design.bands[selected].offset,scale:geometryScale()};}
 else {history.pop();return;}
 svg.setPointerCapture(e.pointerId);panel();
});
svg.addEventListener('pointermove',e=>{
 if(!drag)return;const point=coords(e);
 if(drag.type==='anchor'){
  const b=design.bands[drag.band],count=compiled.counts[drag.band];
  b.offset=clamp(drag.initial-(point.y-drag.y)/drag.scale,-b.radius+.000025,.02);
  const start=design.bands.slice(0,drag.band).reduce((s,b)=>s+b.length,0);
  b.anchor=clamp((((point.x-60)/880*width()-start)/b.length-.5/count)/(1-1/count),0,1);
 }
 edited();
});
const endDrag=()=>{if(drag){drag=null;clearTimeout(editTimer);refresh();}};svg.addEventListener('pointerup',endDrag);svg.addEventListener('pointercancel',endDrag);
svg.addEventListener('keydown',e=>{
 const anchor=e.target.closest('.anchor');if(anchor&&['ArrowUp','ArrowDown','ArrowLeft','ArrowRight'].includes(e.key)){
 e.preventDefault();checkpoint();const b=design.bands[+anchor.dataset.band];
 if(e.key==='ArrowUp'||e.key==='ArrowDown')b.offset=clamp(b.offset+(e.key==='ArrowUp'?1:-1)*.0001,-b.radius+.000025,.02);
 else b.anchor=clamp(b.anchor+(e.key==='ArrowRight'?.03:-.03),0,1);refresh();svg.querySelector(`.anchor[data-band="${anchor.dataset.band}"]`)?.focus();
 }

});
$('bands').onclick=e=>{const b=e.target.closest('button');if(b){selected=+b.dataset.band;panel();draw();}};
for(const k of ['mass','elasticity','support','pressure','wander']){
 $(k).addEventListener('pointerdown',checkpoint);$(k).addEventListener('keydown',e=>{if(e.key.startsWith('Arrow'))checkpoint();});
 $(k).oninput=()=>{if(['pressure','wander'].includes(k))system[k]=+$(k).value;else design.bands[selected][k]=+$(k).value;refresh();};
}
for(const [id,key] of [['body-elasticity','elasticityScale'],['elasticity-variation','elasticityVariation']]){
 $(id).onpointerdown=checkpoint;$(id).onkeydown=e=>{if(e.key.startsWith('Arrow'))checkpoint();};$(id).oninput=()=>{system[key]=+$(id).value;refresh();};
}
function choosePreset(value){checkpoint();system=value==='reference'?referenceDesign():value==='below-mouth'?initialDesign():systemDesign(value);stage=0;systemState=state=null;selected=system.version===4?4:compile(system).cordIndex;modeIndex=0;refresh(false);if(node)node.port.postMessage({type:'reset',compiled:whole});benchmark();}
$('preset').onchange=()=>choosePreset($('preset').value);
$('material').onchange=()=>{checkpoint();design.bands[selected].material=$('material').value;refresh();};
$('monitor').oninput=()=>{if(node)node.port.postMessage({type:'level',value:+$('monitor').value});panel();};
$('links').onchange=draw;$('mode-view').onchange=draw;$('mode-select').onchange=()=>{modeIndex=+$('mode-select').value;previewShape=modeShape(compiled.sheet,modeIndex);$('mode-view').checked=true;draw();};
$('undo').onclick=()=>{if(history.length){system=JSON.parse(history.pop());systemState=state=null;refresh();}};
$('reset').onclick=()=>choosePreset(system.version===4?system.arrangement:system.version===3?'below-mouth':'reference');
function download(data,type,name){const url=URL.createObjectURL(new Blob([data],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);}
$('save').onclick=()=>download(JSON.stringify(system,null,2),'application/json','ducklips-system.json');
$('load').onchange=async e=>{try{const file=e.target.files[0];if(!file)return;if(file.size>100000)throw Error('Design file is too large.');const value=validateSystem(JSON.parse(await file.text()));checkpoint();system=value;stage=0;selected=0;systemState=state=null;refresh();benchmark();notice('System loaded.');}catch(e){notice(e.message);}finally{$('load').value='';}};
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
 node.port.postMessage({type:'configure',compiled:whole});node.port.postMessage({type:'level',value:+$('monitor').value});
 node.port.onmessage=({data:d})=>{if(d.type==='state'){systemState=d;state=d.tubes?.[stage]??d;$('elasticity-live').textContent=state.elasticityFactors?`Selected material: ${((state.elasticityFactors[selected]-1)*100).toFixed(2)}% from its set elasticity`:'';wave=d.scope;waveIndex=d.scopeIndex;$('pressure-reading').textContent=Math.round(d.drive)+' Pa inlet';if(d.fault){running=false;node.port.postMessage({type:'run',value:false});$('start').textContent='Restart pressure';$('status').textContent='Model stopped';notice(d.fault);}}
 if(d.type==='recording'){addTake(d.samples,d.rate,'Live tube',captureContext);captureContext=null;$('record').disabled=false;$('record').textContent='Record 2 seconds';}
 };
 await context.resume();
}
$('start').onclick=async()=>{try{await audio();if(systemState?.fault){node.port.postMessage({type:'reset',compiled:whole});state=null;}running=!running;node.port.postMessage({type:'run',value:running});$('start').textContent=running?'Stop pressure':'Start pressure';$('status').textContent=running?'Running · outlet connected':'Pressure off · releasing';}catch(e){notice('Audio could not start: '+e.message);}};
$('record').onclick=async()=>{try{await audio();if(!running){running=true;node.port.postMessage({type:'run',value:true});$('start').textContent='Stop pressure';$('status').textContent='Running · recording outlet';}$('record').disabled=true;$('record').textContent='Recording…';captureContext={started:performance.now(),date:new Date().toISOString(),design:JSON.parse(JSON.stringify(system)),initialObservation:systemState,monitorGain:+$('monitor').value,changes:[]};node.port.postMessage({type:'capture',seconds:2});}catch(e){notice(e.message);$('record').disabled=false;}};
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
let benchmarkWorker;
function benchmark(){benchmarkWorker?.terminate();$('speed').textContent='Benchmarking…';const worker=benchmarkWorker=new Worker('./render-worker.mjs',{type:'module'});worker.onmessage=({data:d})=>{if(d.error){notice(d.error);return;}$('speed').textContent=d.realtime.toFixed(1)+'× core realtime';if(d.realtime<10)notice('Core is below the 10× target on this browser.');worker.terminate();};worker.postMessage({compiled:whole,seconds:3,rate:24000});}
benchmark();
