import http from 'node:http';
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {spawn} from 'node:child_process';
import {createHash} from 'node:crypto';
import {networkInterfaces} from 'node:os';
const root=path.dirname(fileURLToPath(import.meta.url)),repo=path.resolve(root,'../..'),out=path.join(repo,'out/flat-tube-web');fs.mkdirSync(out,{recursive:true});
const types={'.html':'text/html','.css':'text/css','.mjs':'text/javascript','.js':'text/javascript','.json':'application/json','.png':'image/png','.wav':'audio/wav','.npz':'application/octet-stream','.md':'text/plain'};let busy=false;
const server=http.createServer(async(req,res)=>{
 let pathname;try{pathname=decodeURIComponent(new URL(req.url,'http://localhost').pathname);}catch{res.writeHead(400);res.end('Invalid URL');return;}
 const json=(status,data)=>{res.writeHead(status,{'Content-Type':'application/json'});res.end(JSON.stringify(data));};
 if(pathname==='/api/enhanced'&&req.method==='POST'){
  if(busy)return json(409,{error:'An enhanced field is already computing. Retry when it finishes.'});
  if(req.headers.origin&&req.headers.origin!==`http://${req.headers.host}`)return json(403,{error:'Local page origin required.'});
  busy=true;try{
   const chunks=[];let size=0;for await(const chunk of req){size+=chunk.length;if(size>2_000_000)throw Error('Recording exceeds the 2 MB limit.');chunks.push(chunk);}let data=Buffer.concat(chunks),provenance=null;if(req.headers['content-type']==='application/json'){const parsed=JSON.parse(data.toString());if(typeof parsed.audio!=='string')throw Error('Expected audio.');data=Buffer.from(parsed.audio,'base64');provenance=parsed.provenance;}
   if(data.length<44||data.toString('ascii',0,4)!=='RIFF'||data.toString('ascii',8,12)!=='WAVE')throw Error('Expected WAV audio.');
   const id=new Date().toISOString().replaceAll(':','-')+'-'+Math.random().toString(16).slice(2,8),dir=path.join(out,id);fs.mkdirSync(dir);const wav=path.join(dir,'outlet.wav');fs.writeFileSync(wav,data);
   const sourceFiles=['engine.mjs','sheet.mjs','audio-renderer.mjs','worklet.mjs','lan-audio.mjs','audio-worker.mjs','app.mjs'];fs.mkdirSync(path.join(dir,'source'));const hashes={};for(const file of sourceFiles){const bytes=fs.readFileSync(path.join(root,file));fs.writeFileSync(path.join(dir,'source',file),bytes);hashes[file]=createHash('sha256').update(bytes).digest('hex');}
   fs.writeFileSync(path.join(dir,'take.json'),JSON.stringify({audio_sha256:createHash('sha256').update(data).digest('hex'),source_sha256:hashes,provenance},null,2));
   const python=process.env.DUCKLIPS_PYTHON||path.join(process.env.HOME,'miniforge3/bin/python');
   const child=spawn(python,[path.join(repo,'tools/enhanced_specgram/render.py'),wav,'--channel','0','--max-hz','8000','--hop','256','--out',path.join(dir,'enhanced')],{env:{...process.env,OPENBLAS_NUM_THREADS:'1'}});
   let log='';child.stdout.on('data',b=>{log+=b;});child.stderr.on('data',b=>{log+=b;});
   await new Promise((resolve,reject)=>{child.once('error',reject);child.once('close',code=>code===0?resolve():reject(Error('Enhanced renderer failed; see local render log.')));});
   fs.writeFileSync(path.join(dir,'render.log'),log);fs.writeFileSync(path.join(dir,'index.html'),`<!doctype html><meta charset="utf-8"><title>Tube · full enhanced field</title><style>body{background:#08090a;color:#ddd;font:14px system-ui;margin:24px}img{display:block;max-width:none;image-rendering:pixelated}p{max-width:800px}a{color:#d8ed96}</style><h1>Complete enhanced field</h1><p>Every computed cell, original magma. Scroll to inspect the full image. No ordinary STFT companion, crop, or image reduction.</p><audio controls src="outlet.wav"></audio><p><a href="enhanced.json">Provenance</a> · <a href="enhanced.npz">Field data</a></p><img src="enhanced.png" alt="Full enhanced depiction of this exact outlet recording">`);
   json(200,{url:`/renders/${id}/index.html`});
  }catch(e){json(400,{error:e.message});}finally{busy=false;}return;
 }
 if(req.method!=='GET'&&req.method!=='HEAD'){res.writeHead(405);res.end();return;}
 let base=root,relative=pathname==='/'?'index.html':pathname.slice(1);if(pathname.startsWith('/renders/')){base=out;relative=pathname.slice(9);}
 const file=path.resolve(base,relative);if(!file.startsWith(base+path.sep)){res.writeHead(403);res.end();return;}
 try{const stat=fs.statSync(file);if(!stat.isFile())throw Error();res.writeHead(200,{'Content-Type':types[path.extname(file)]||'application/octet-stream','Cache-Control':'no-store','Content-Length':stat.size});if(req.method==='HEAD')res.end();else fs.createReadStream(file).pipe(res);}catch{res.writeHead(404);res.end('Not found');}
});
const host=process.env.HOST||'0.0.0.0';
server.listen(Number(process.env.PORT||8765),host,()=>{
 const port=server.address().port;
 console.log(`Ducklips tube workshop: http://127.0.0.1:${port}/`);
 if(host==='0.0.0.0')for(const list of Object.values(networkInterfaces()))for(const address of list||[])if(address.family==='IPv4'&&!address.internal)console.log(`LAN: http://${address.address}:${port}/`);
});
