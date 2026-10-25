"""Sliders for the owner: a small local web page that renders the engine
with the settings the sliders hold.

The page (served at http://localhost:PORT) shows one slider per setting of
the pinch tube and its neighbours (posture, breath noise, the tract's
friction), a choice of utterance (held ah / ee / oo, the banana sentence),
and plays each render. Renders run here (speak, FOLDTUBE=1, PARAMS = the
base tree with the sliders' values written in); presets are saved under
tract/data/sliders/NAME.params with the values as JSON beside them.

usage: slider_server.py BASE.params [PORT]   (bound to localhost only)"""
import http.server
import json
import os
import re
import subprocess
import sys
import tempfile
import threading
import time

ROOT = os.path.expanduser("~/ducklips")
SPEAK = os.path.join(ROOT, "tract", "build", "speak")
SCORES = os.path.join(ROOT, "tract", "data", "slider_scores")
PRESETS = os.path.join(ROOT, "tract", "data", "sliders")
UTTERANCES = {"ah": "ah.score", "ee": "ee.score", "oo": "oo.score", "banana": "banana.score"}

# (block.key, label, low, high, log scale, group)
SLIDERS = [
    ("pinch.support", "support (wall to rest)", 1e4, 3e5, True, "Pinch walls"),
    ("pinch.elasticity", "elasticity (wall to wall)", 1e3, 1e5, True, "Pinch walls"),
    ("pinch.surface_density", "mass (surface density)", 0.03, 0.5, True, "Pinch walls"),
    ("pinch.loss", "damping", 0.02, 0.5, True, "Pinch walls"),
    ("pinch.contact_stiffness", "contact stiffness", 3e4, 3e6, True, "Pinch walls"),
    ("pinch.contact_loss", "contact damping", 0.1, 2.0, False, "Pinch walls"),
    ("pinch.natural_pitch", "tension anchor (natural pitch)", 60, 250, False, "Pinch walls"),
    ("pinch.width", "width across (cm)", 0.6, 2.0, False, "Geometry"),
    ("pinch.section_length", "section length (cm)", 0.1, 1.0, False, "Geometry"),
    ("pinch.outer_gap", "outer gap, sections 1 & 6 (cm)", 0.2, 1.5, False, "Geometry"),
    ("pinch.inner_gap", "inner gap, sections 2 & 5 (cm)", 0.05, 0.8, False, "Geometry"),
    ("posture.voicing_lower", "pinch rest, section 3 (cm, each side)", -0.01, 0.03, False, "Geometry"),
    ("posture.voicing_upper", "pinch rest, section 4 (cm, each side)", -0.01, 0.03, False, "Geometry"),
    ("posture.voicing_chink", "leak beside the pinch (cm2)", 0.0, 0.02, False, "Geometry"),
    ("pinch.wander", "elasticity wander", 0.0, 0.3, False, "Wander"),
    ("pinch.wander_corner", "wander speed (Hz)", 0.5, 30.0, True, "Wander"),
    ("pinch.wander_tamping", "wander tamped by support", 0.0, 3.0, False, "Wander"),
    ("pinch.jet_turbulence", "jet eddies on the walls (share of jet pressure)", 0.0, 1.0, False, "Jet"),
    ("pinch.turbulence_corner", "eddy speed (Hz)", 100.0, 5000.0, True, "Jet"),
    ("pinch.string_modes", "string modes (0 = piston walls)", 0.0, 12.0, False, "Harp strings"),
    ("pinch.string_tension", "string tension across the slit (dyn)", 1000.0, 100000.0, True, "Harp strings"),
    ("pinch.string_loss_growth", "higher modes die faster (loss n^x)", 0.0, 2.0, False, "Harp strings"),
    ("pinch.string_taper", "rest gap opening toward one end (cm)", 0.0, 0.05, False, "Harp strings"),
]
for section in range(1, 7):
    for piece in ("mass", "elasticity", "support"):
        SLIDERS.append(("pinch.s%d_%s" % (section, piece), "section %d %s" % (section, piece), 0.02, 10.0, True,
                        "Section %d" % section))
SLIDERS += [
    ("glottal_noise.aspiration", "breath noise", 0.0, 40.0, False, "Air"),
    ("glottal_noise.fold_turbulence", "turbulence noise in the sound", 0.0, 10.0, False, "Air"),
    ("tract.viscous_scale", "tract friction", 0.05, 1.0, True, "Air"),
    ("walls.loss_scale", "tract wall loss", 0.1, 2.0, True, "Air"),
]

LOCK = threading.Lock()


def read_value(text, leaf):
    block, key = leaf.split(".")
    inside = False
    found = None
    for line in text.split("\n"):
        if re.match(r"^%s \{" % re.escape(block), line):
            inside = True
        elif inside and line.startswith("}"):
            inside = False
        elif inside and line.strip().split()[:1] == [key]:
            found = float(line.split()[1])
    return found


def with_values(text, values):
    """The tree with each block.key set (the last block of that name holds
    it; a missing key is added to it; a missing block is appended)."""
    lines = text.rstrip("\n").split("\n")
    for leaf, value in values.items():
        block, key = leaf.split(".")
        starts = [i for i, line in enumerate(lines) if re.match(r"^%s \{" % re.escape(block), line)]
        if not starts:
            lines += ["%s {" % block, "  %s %.10g" % (key, value), "}"]
            continue
        start = starts[-1]
        end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("}"))
        for i in range(start + 1, end):
            if lines[i].strip().split()[:1] == [key]:
                lines[i] = "  %s %.10g" % (key, value)
                break
        else:
            lines.insert(end, "  %s %.10g" % (key, value))
    return "\n".join(lines) + "\n"


def defaults(base_text):
    out = {}
    for leaf, label, low, high, log, group in SLIDERS:
        value = read_value(base_text, leaf)
        out[leaf] = value if value is not None else low
    return out


PAGE = r"""<!doctype html><html><head><meta charset="utf-8"><title>Pinch Sliders</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
:root{--bg:#f6f4ef;--card:#fff;--ink:#222;--dim:#777;--accent:#b5523b;--line:#ddd}
@media (prefers-color-scheme: dark){:root{--bg:#1b1a18;--card:#262421;--ink:#eee;--dim:#999;--accent:#e0805f;--line:#3a3733}}
body{background:var(--bg);color:var(--ink);font:14px/1.35 -apple-system,system-ui,sans-serif;margin:0;padding:12px 16px}
h1{font-size:17px;margin:4px 0 10px}
.bar{position:sticky;top:0;background:var(--bg);padding:8px 0;z-index:2;display:flex;flex-wrap:wrap;gap:8px;align-items:center;border-bottom:1px solid var(--line)}
button,select,input[type=text]{font:inherit;padding:5px 9px;border-radius:6px;border:1px solid var(--line);background:var(--card);color:var(--ink)}
button.go{background:var(--accent);color:#fff;border-color:var(--accent)}
.groups{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:10px;margin-top:10px}
.group{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px 10px}
.group h2{font-size:13px;margin:0 0 6px;color:var(--accent)}
.row{display:grid;grid-template-columns:1fr 70px;gap:2px 8px;align-items:center;margin:4px 0}
.row label{font-size:12px;color:var(--dim)}
.row input[type=range]{grid-column:1/3;width:100%}
.val{font-variant-numeric:tabular-nums;text-align:right;font-size:12px}
#takes{display:flex;flex-direction:column;gap:4px;margin-top:8px}
.take{display:flex;gap:8px;align-items:center;font-size:12px}
.take audio{height:28px}
#status{color:var(--dim);font-size:12px}
</style></head><body>
<h1>Pinch tube sliders</h1>
<div class="bar">
 <select id="utt"><option value="ah">held ah</option><option value="ee">held ee</option><option value="oo">held oo</option><option value="banana">banana sentence</option></select>
 <button class="go" id="render">Render</button>
 <label><input type="checkbox" id="auto" checked> render on release</label>
 <button id="reset">Reset</button>
 <input type="text" id="name" placeholder="preset name" size="12"><button id="save">Save preset</button>
 <span id="status"></span>
</div>
<div id="takes"></div>
<div class="groups" id="groups"></div>
<script>
let specs=[], values={}, base={}, count=0;
const toSlider=(s,v)=>s.log?1000*Math.log(v/s.low)/Math.log(s.high/s.low):1000*(v-s.low)/(s.high-s.low);
const fromSlider=(s,x)=>s.log?s.low*Math.pow(s.high/s.low,x/1000):s.low+(s.high-s.low)*x/1000;
const fmt=v=>Math.abs(v)>=1000?v.toExponential(2):(Math.abs(v)<0.01&&v!==0?v.toExponential(2):(+v.toPrecision(3)).toString());
async function load(){
 const r=await fetch('/spec'); const j=await r.json(); specs=j.sliders; base=j.values; values=Object.assign({},base); build();
}
function build(){
 const g=document.getElementById('groups'); g.innerHTML=''; const groups={};
 for(const s of specs){
  if(!groups[s.group]){const d=document.createElement('div');d.className='group';d.innerHTML='<h2>'+s.group+'</h2>';g.appendChild(d);groups[s.group]=d;}
  const row=document.createElement('div');row.className='row';
  row.innerHTML='<label>'+s.label+'</label><span class="val"></span><input type="range" min="0" max="1000" step="1">';
  const input=row.querySelector('input'), val=row.querySelector('.val');
  input.value=toSlider(s,values[s.leaf]); val.textContent=fmt(values[s.leaf]);
  input.oninput=()=>{values[s.leaf]=fromSlider(s,+input.value);val.textContent=fmt(values[s.leaf]);};
  input.onchange=()=>{if(document.getElementById('auto').checked)render();};
  groups[s.group].appendChild(row);
 }
}
async function render(){
 const utt=document.getElementById('utt').value, st=document.getElementById('status');
 st.textContent='rendering…'; const t0=performance.now();
 const r=await fetch('/render',{method:'POST',body:JSON.stringify({utterance:utt,values:values})});
 if(!r.ok){st.textContent='render failed: '+await r.text();return;}
 const blob=await r.blob(); const url=URL.createObjectURL(blob); count++;
 const take=document.createElement('div');take.className='take';
 take.innerHTML='<b>#'+count+'</b> '+utt+' <audio controls src="'+url+'"></audio>';
 const takes=document.getElementById('takes'); takes.prepend(take);
 while(takes.children.length>6)takes.removeChild(takes.lastChild);
 take.querySelector('audio').play();
 st.textContent='take #'+count+' in '+((performance.now()-t0)/1000).toFixed(1)+' s';
}
document.getElementById('render').onclick=render;
document.getElementById('reset').onclick=()=>{values=Object.assign({},base);build();};
document.getElementById('save').onclick=async()=>{
 const name=document.getElementById('name').value.trim(); const st=document.getElementById('status');
 if(!/^[A-Za-z0-9_-]+$/.test(name)){st.textContent='name: letters, digits, - or _';return;}
 const r=await fetch('/save',{method:'POST',body:JSON.stringify({name:name,values:values})});
 st.textContent=r.ok?'saved '+name:'save failed';
};
load();
</script></body></html>"""


class Handler(http.server.BaseHTTPRequestHandler):
    base_text = ""

    def send(self, code, body, kind):
        self.send_response(code)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/":
            self.send(200, PAGE.encode(), "text/html; charset=utf-8")
        elif self.path == "/spec":
            sliders = [dict(leaf=leaf, label=label, low=low, high=high, log=log, group=group)
                       for leaf, label, low, high, log, group in SLIDERS]
            self.send(200, json.dumps(dict(sliders=sliders, values=defaults(self.base_text))).encode(),
                      "application/json")
        else:
            self.send(404, b"not found", "text/plain")

    def do_POST(self):
        length = int(self.headers.get("Content-Length", "0"))
        request = json.loads(self.rfile.read(length) or b"{}")
        known = {leaf for leaf, _, _, _, _, _ in SLIDERS}
        values = {leaf: float(value) for leaf, value in request.get("values", {}).items() if leaf in known}
        if self.path == "/render":
            utterance = request.get("utterance", "ah")
            if utterance not in UTTERANCES:
                self.send(400, b"unknown utterance", "text/plain")
                return
            with LOCK:
                folder = tempfile.mkdtemp(prefix="sliders_")
                params = os.path.join(folder, "s.params")
                wave = os.path.join(folder, "s.wav")
                with open(params, "w") as handle:
                    handle.write(with_values(self.base_text, values))
                env = dict(os.environ, FOLDTUBE="1", PARAMS=params)
                done = subprocess.run([SPEAK, os.path.join(SCORES, UTTERANCES[utterance]), wave, "73021"], env=env,
                                      capture_output=True, text=True, timeout=120)
                if done.returncode != 0 or not os.path.exists(wave):
                    self.send(500, (done.stderr[-400:] or "speak failed").encode(), "text/plain")
                    return
                body = open(wave, "rb").read()
                for path in (params, wave):
                    os.remove(path)
                os.rmdir(folder)
            log = os.path.join(PRESETS, "render_log.jsonl")
            with open(log, "a") as handle:
                handle.write(json.dumps(dict(time=time.time(), utterance=utterance, values=values)) + "\n")
            self.send(200, body, "audio/wav")
        elif self.path == "/save":
            name = str(request.get("name", ""))
            if not re.match(r"^[A-Za-z0-9_-]+$", name):
                self.send(400, b"bad name", "text/plain")
                return
            with open(os.path.join(PRESETS, name + ".params"), "w") as handle:
                handle.write("# Slider preset %s (%s)\n" % (name, time.strftime("%Y-%m-%d %H:%M")) +
                             with_values(self.base_text, values))
            with open(os.path.join(PRESETS, name + ".json"), "w") as handle:
                json.dump(values, handle, indent=1)
            self.send(200, b"ok", "text/plain")
        else:
            self.send(404, b"not found", "text/plain")

    def log_message(self, format, *args):
        return


def canonical_tree(base_path):
    """The engine's whole tree for the base (every setting with its value,
    the defaults included): speak's PARAMS_OUT."""
    folder = tempfile.mkdtemp(prefix="sliders_")
    out = os.path.join(folder, "full.params")
    env = dict(os.environ, FOLDTUBE="1", PARAMS=base_path, PARAMS_OUT=out)
    subprocess.run([SPEAK, os.path.join(SCORES, UTTERANCES["ah"]), os.path.join(folder, "x.wav"), "73021"], env=env,
                   capture_output=True, timeout=120)
    text = open(out).read()
    for name in os.listdir(folder):
        os.remove(os.path.join(folder, name))
    os.rmdir(folder)
    return text


def main():
    base_path = sys.argv[1]
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 8765
    os.makedirs(PRESETS, exist_ok=True)
    Handler.base_text = canonical_tree(base_path)
    server = http.server.ThreadingHTTPServer(("127.0.0.1", port), Handler)
    print("sliders at http://localhost:%d (base %s)" % (port, base_path), flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
