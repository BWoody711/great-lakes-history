"use strict";
(function(){

/* ---------------------------------------------------------------------
 * Data loading. Both datasets are fetched at runtime rather than embedded,
 * so the page stays small and the JSON stays a single source of truth
 * shared with the standalone viewers/*.html archives.
 * ------------------------------------------------------------------- */
const DATA_FILES = {
  optimal: 'data/paleo_optimal.json',
  max:     'data/paleo_max.json',
  algonquin: 'data/algonquin.json',
  rebound: 'data/rebound.json',
};

Promise.all(Object.entries(DATA_FILES).map(([k,u]) =>
  fetch(u).then(r => { if(!r.ok) throw new Error(u+': '+r.status); return r.json(); }).then(j => [k,j])
)).then(entries => {
  const D = Object.fromEntries(entries);
  init({optimal: D.optimal, max: D.max}, D.algonquin, D.rebound);
}).catch(err => {
  document.getElementById('loading').textContent = 'failed to load data: '+err.message;
});

function init(LK_DS, ALG, RB_D){
"use strict";

/* ---- lake / palaeogeography state ---- */
let LK_MARGIN='optimal';
let LK_D=LK_DS[LK_MARGIN];
let LK_F=LK_D.frames.slice().sort((a,b)=>b.t-a.t);
let MODERN=LK_F[LK_F.length-1];

/* ---- rebound state ---- */
const RB_F=RB_D.frames.slice().sort((a,b)=>b.t-a.t);
const RB_LON=RB_D.lon1d, RB_LAT=RB_D.lat1d, RB_BG=RB_D.bg;
let RB_MODE='none', RB_MODEL='ice7';

/* ---- shared map viewport ---- */
const VW=-94, VE=-71, VN=51.5, VS=39;
const cv=document.getElementById('map'), ctx=cv.getContext('2d');
const CW=cv.width, CH=cv.height;
const px = lo => (lo-VW)/(VE-VW)*CW, py = la => (VN-la)/(VN-VS)*CH;

document.getElementById('loading').hidden = true;

/* ================= lake-overlay drawing (Path2D, cached) ================= */
function pathOf(rings){
  const p=new Path2D();
  for(const poly of rings){
    for(const ring of poly){
      p.moveTo(px(ring[0][0]), py(ring[0][1]));
      for(let i=1;i<ring.length;i++) p.lineTo(px(ring[i][0]), py(ring[i][1]));
      p.closePath();
    }
  }
  return p;
}
const pcache=new Map();
function cached(key, rings){
  key=LK_MARGIN+key;
  let p=pcache.get(key); if(!p){ p=pathOf(rings); pcache.set(key,p);} return p;
}

/* ================= rebound raster background ================= */
function giaAt(grid, lo, la){
  let i=Math.floor(lo-RB_LON[0]), j=Math.floor(la-RB_LAT[0]);
  if(i<0||j<0||i>=RB_LON.length-1||j>=RB_LAT.length-1) return null;
  const fx=lo-RB_LON[0]-i, fy=la-RB_LAT[0]-j;
  const a=grid[j][i], b=grid[j][i+1], c=grid[j+1][i], d=grid[j+1][i+1];
  if(a==null||b==null||c==null||d==null) return null;
  return (a*(1-fx)+b*fx)*(1-fy)+(c*(1-fx)+d*fx)*fy;
}
function gridFor(t,model){
  const src = model==='ice7' ? RB_D.ice7 : RB_D.ice6;
  if(t===0) return null;
  return src[String(t)] || src[String(nearestT(t))];
}
function upliftColour(v,max){
  const x=Math.max(-1,Math.min(1,v/max));
  if(x>=0){ const u=Math.pow(x,.75);
    return [Math.round(250-140*u), Math.round(246-120*u), Math.round(225-40*u)]; }
  const u=Math.pow(-x,.75);
  return [Math.round(250-30*u), Math.round(246-120*u), Math.round(225-150*u)];
}
function topoColour(z){
  if(z<0){ const u=Math.min(1,-z/350);
    return [Math.round(120-70*u),Math.round(165-70*u),Math.round(200-45*u)]; }
  const stops=[[0,[186,205,170]],[150,[214,214,168]],[300,[203,186,138]],
               [600,[176,150,112]],[1200,[150,132,116]],[2500,[235,235,235]]];
  for(let i=0;i<stops.length-1;i++){
    const [z0,c0]=stops[i],[z1,c1]=stops[i+1];
    if(z<=z1){const u=(z-z0)/(z1-z0);
      return [0,1,2].map(k=>Math.round(c0[k]+(c1[k]-c0[k])*u));}
  }
  return [235,235,235];
}
let img=null;
function drawRaster(t){
  const g=gridFor(t,RB_MODEL), g6=gridFor(t,'ice6'), g7=gridFor(t,'ice7');
  if(!img) img=ctx.createImageData(CW,CH);
  const px0=img.data;
  const bw=RB_BG.w, bh=RB_BG.h;
  const [W,S,E,N]=RB_BG.bbox;
  let smax = RB_MODE==='diff' ? 60 : 500;
  if(RB_MODE==='uplift'){ smax=0; if(g) for(const row of g) for(const v of row) if(v!=null) smax=Math.max(smax,Math.abs(v)); smax=Math.max(50,Math.ceil(smax/50)*50); }
  for(let y=0;y<CH;y++){
    const la = VN-(y+0.5)/CH*(VN-VS);
    const by = Math.floor((N-la)/(N-S)*bh);
    for(let x=0;x<CW;x++){
      const lo = VW+(x+0.5)/CW*(VE-VW);
      const bx = Math.floor((lo-W)/(E-W)*bw);
      const k=(y*CW+x)*4;
      let r,gg,b;
      if(bx<0||by<0||bx>=bw||by>=bh){ r=gg=b=20; }
      else {
        const z0=RB_BG.z[by*bw+bx];
        if(RB_MODE==='topo'){
          const u = t===0?0:(giaAt(g,lo,la)||0);
          [r,gg,b]=topoColour(z0-u);
        } else if(RB_MODE==='diff'){
          const a=t===0?0:(giaAt(g7,lo,la)||0), c=t===0?0:(giaAt(g6,lo,la)||0);
          [r,gg,b]=upliftColour(a-c,smax);
        } else {
          const u=t===0?0:(giaAt(g,lo,la)||0);
          [r,gg,b]=upliftColour(u,smax);
          if(z0<0){ r=Math.round(r*.72+30); gg=Math.round(gg*.72+40); b=Math.round(b*.72+60); }
        }
      }
      px0[k]=r; px0[k+1]=gg; px0[k+2]=b; px0[k+3]=255;
    }
  }
  ctx.putImageData(img,0,0);
  if(t!==0){
    ctx.lineWidth=1; ctx.strokeStyle='rgba(30,40,50,.35)';
    const step = RB_MODE==='diff'?20:100;
    for(let lvl=-600;lvl<=800;lvl+=step){
      ctx.beginPath();
      for(let x=0;x<CW;x+=3){
        const lo=VW+x/CW*(VE-VW);
        let prev=null;
        for(let y=0;y<CH;y+=3){
          const la=VN-y/CH*(VN-VS);
          const v = RB_MODE==='diff' ? (giaAt(g7,lo,la)-giaAt(g6,lo,la)) : giaAt(g,lo,la);
          if(v==null){prev=null;continue;}
          if(prev!==null && ((prev<lvl)!==(v<lvl))) { ctx.moveTo(x,y); ctx.lineTo(x+3,y); }
          prev=v;
        }
      }
      ctx.stroke();
    }
  }
  drawRamp(smax);
}
function drawRamp(smax){
  const rc=document.getElementById('ramp'), rx=rc.getContext('2d');
  const id=rx.createImageData(150,11);
  for(let x=0;x<150;x++){
    let c;
    if(RB_MODE==='topo') c=topoColour(-200+ x/149*1400);
    else c=upliftColour(-smax + x/149*2*smax, smax);
    for(let y=0;y<11;y++){const k=(y*150+x)*4; id.data[k]=c[0];id.data[k+1]=c[1];id.data[k+2]=c[2];id.data[k+3]=255;}
  }
  rx.putImageData(id,0,0);
  document.getElementById('smin').textContent = RB_MODE==='topo' ? '−200 m' : (-smax)+' m';
  document.getElementById('smax').textContent = RB_MODE==='topo' ? '1200 m' : '+'+smax+' m';
}

/* ================= combined map draw ================= */
function drawMap(){
  const t=nearestT(cur);
  const lkfr=nearest(cur);

  if(RB_MODE==='none'){
    ctx.fillStyle='#b4b2a1'; ctx.fillRect(0,0,CW,CH);
    ctx.strokeStyle='rgba(0,0,0,.07)'; ctx.lineWidth=1;
    for(let lo=-92;lo<=-74;lo+=2){ctx.beginPath();ctx.moveTo(px(lo),0);ctx.lineTo(px(lo),CH);ctx.stroke();}
    for(let la=41;la<=50;la++){ctx.beginPath();ctx.moveTo(0,py(la));ctx.lineTo(CW,py(la));ctx.stroke();}
    document.getElementById('scaleWrap').style.display='none';
  } else {
    drawRaster(t);
    document.getElementById('scaleWrap').style.display='flex';
  }

  if(document.getElementById('lyTransect').checked){
    const Tr=RB_D.transect;
    ctx.strokeStyle='#b1372a'; ctx.lineWidth=2.5; ctx.setLineDash([9,5]);
    ctx.beginPath(); ctx.moveTo(px(Tr.lon[0]),py(Tr.lat[0]));
    ctx.lineTo(px(Tr.lon[Tr.lon.length-1]),py(Tr.lat[Tr.lat.length-1])); ctx.stroke();
    ctx.setLineDash([]);
    for(const m of Tr.marks){
      ctx.beginPath(); ctx.arc(px(m.lon),py(m.lat),4,0,6.284);
      ctx.fillStyle='#fff'; ctx.fill(); ctx.strokeStyle='#b1372a'; ctx.lineWidth=2; ctx.stroke();
    }
  }

  if(document.getElementById('lyModern').checked){
    ctx.save(); ctx.setLineDash([7,6]); ctx.strokeStyle='rgba(60,64,50,.55)'; ctx.lineWidth=1.4;
    for(const lk of MODERN.lakes) if(lk.a>3000) ctx.stroke(cached('m'+lk.lv+lk.a, lk.r));
    ctx.restore();
  }
  if(document.getElementById('lySea').checked && lkfr.sea.length){
    ctx.fillStyle='#3a8b7e'; ctx.fill(cached('s'+lkfr.t, lkfr.sea),'evenodd');
  }
  if(document.getElementById('lyLakes').checked){
    for(const lk of lkfr.lakes){
      ctx.fillStyle = lk.a>4000 ? '#22537d' : '#5f93ba';
      ctx.fill(cached('l'+lkfr.t+'_'+lk.lv+'_'+lk.a, lk.r),'evenodd');
    }
  }
  if(document.getElementById('lyIce').checked && lkfr.ice.length){
    const p=cached('i'+lkfr.t, lkfr.ice);
    ctx.fillStyle='#e9f1f4'; ctx.fill(p,'evenodd');
    ctx.strokeStyle='#93b4c3'; ctx.lineWidth=2; ctx.stroke(p);
  }
  if(document.getElementById('lySpill').checked && document.getElementById('lyLakes').checked){
    ctx.lineWidth=2.5;
    for(const lk of lkfr.lakes){
      if(!lk.sp || lk.a<8000) continue;
      const x=px(lk.sp[0]), y=py(lk.sp[1]);
      ctx.beginPath(); ctx.arc(x,y,6,0,6.284);
      ctx.fillStyle='#fff'; ctx.fill(); ctx.strokeStyle='#a24334'; ctx.stroke();
    }
  }
  if(document.getElementById('lyAlg').checked){
    const a=(ALG.frames[LK_MARGIN]||[]).find(x=>x.t===lkfr.t);
    if(a && a.r.length){
      ctx.save();
      ctx.setLineDash([7,5]); ctx.lineWidth=2.6; ctx.strokeStyle='#e0512f';
      ctx.stroke(cached('a'+lkfr.t, a.r));
      ctx.restore();
      ctx.font='italic 19px Georgia, serif'; ctx.textAlign='left';
      ctx.lineWidth=4; ctx.strokeStyle='#2a1410';
      const lab='Lake Algonquin, '+ALG.strandline+' m strandline';
      ctx.strokeText(lab, 26, CH-26); ctx.fillStyle='#f0a08c';
      ctx.fillText(lab, 26, CH-26);
    }
  }
  if(document.getElementById('lyLabels').checked && document.getElementById('lyLakes').checked){
    ctx.font='italic 22px Georgia, serif'; ctx.textAlign='center';
    for(const lk of lkfr.lakes){
      if(lk.a<9000 || !lk.b.length) continue;
      let cx=0,cy=0,n=0, big=null, ba=0;
      for(const poly of lk.r){ const a=Math.abs(poly[0].length); if(a>ba){ba=a; big=poly[0];} }
      if(!big) continue;
      for(const p of big){cx+=p[0];cy+=p[1];n++;}
      cx/=n; cy/=n;
      const txt = lk.b.join(' + ');
      ctx.lineWidth=4; ctx.strokeStyle='#123146'; ctx.strokeText(txt, px(cx), py(cy));
      ctx.fillStyle='#f2f7fa'; ctx.fillText(txt, px(cx), py(cy));
    }
  }
  return {t, lkfr};
}

/* ---- hydrograph (lake levels by basin) ---- */
const BAS=[{k:'Superior',c:'#7fb6dc'},{k:'Huron',c:'#e0b862'},
           {k:'Erie',c:'#7fbfa8'},{k:'Ontario',c:'#c38fb8'}];
function levelOf(fr,b){
  let best=null;
  for(const lk of fr.lakes) if(lk.b.indexOf(b)>=0 && (!best||lk.a>best.a)) best=lk;
  return best ? best.lv : null;
}
const hy=document.getElementById('hydro');
let hyHead,hyHd;
const H_W=1000,H_H=200,H_PL=46,H_PR=12,H_PT=14,H_PB=26;
let LMIN,LMAX,hx,hv;
function buildHydro(){
  LMIN=1e9; LMAX=-1e9;
  for(const fr of LK_F) for(const b of BAS){const v=levelOf(fr,b.k); if(v!=null){LMIN=Math.min(LMIN,v);LMAX=Math.max(LMAX,v);} }
  LMIN=Math.floor(LMIN/50)*50; LMAX=Math.ceil(LMAX/50)*50;
  hx = t => H_PL + (14500-t)/14500*(H_W-H_PL-H_PR);
  hv = m => H_PT + (LMAX-m)/(LMAX-LMIN)*(H_H-H_PT-H_PB);
  let h='';
  for(let m=LMIN;m<=LMAX;m+=50){
    h+=`<line x1="${H_PL}" y1="${hv(m)}" x2="${H_W-H_PR}" y2="${hv(m)}" stroke="#262d34"/>`;
    h+=`<text x="${H_PL-7}" y="${hv(m)+4}" fill="#71767c" font-size="11" text-anchor="end" font-family="ui-monospace,monospace">${m}</text>`;
  }
  for(const ka of [14,12,10,8,6,4,2,0])
    h+=`<text x="${hx(ka*1000)}" y="${H_H-8}" fill="#71767c" font-size="11" text-anchor="middle" font-family="ui-monospace,monospace">${ka?ka+' ka':'now'}</text>`;
  for(const b of BAS){
    let d='',open=false;
    for(const fr of LK_F){
      const v=levelOf(fr,b.k);
      if(v==null){open=false;continue;}
      d+=(open?'L':'M')+hx(fr.t).toFixed(1)+','+hv(v).toFixed(1)+' '; open=true;
    }
    h+=`<path d="${d}" fill="none" stroke="${b.c}" stroke-width="2" stroke-linejoin="round"/>`;
  }
  let lx=H_PL+8;
  for(const b of BAS){
    h+=`<rect x="${lx}" y="${H_PT+1}" width="13" height="3" fill="${b.c}"/>`;
    h+=`<text x="${lx+18}" y="${H_PT+7}" fill="#9aa0a6" font-size="11.5" font-family="ui-sans-serif,system-ui,sans-serif">${b.k}</text>`;
    lx+=30+b.k.length*7;
  }
  h+=`<line id="hyHead" x1="0" y1="${H_PT-5}" x2="0" y2="${H_H-H_PB}" stroke="#f0ece0" stroke-width="1.5"/>`;
  h+=`<circle id="hyHd" cx="0" cy="${H_PT-5}" r="3.4" fill="#f0ece0"/>`;
  hy.innerHTML=h;
  hyHead=hy.querySelector('#hyHead'); hyHd=hy.querySelector('#hyHd');
}

/* ---- validation table (lakes) ---- */
const REAL={'Superior':[82100,183.5],'Michigan+Huron+Georgian Bay':[117400,176.0],
            'Erie':[25700,174.0],'Ontario':[18960,74.2]};
function buildValid(){
  let h='<tr><th>Basin</th><th style="text-align:right">model</th><th style="text-align:right">actual</th></tr>';
  for(const lk of MODERN.lakes.slice().sort((a,b)=>b.a-a.a)){
    const key=lk.b.join('+');
    if(!(key in REAL)) continue;
    const [ra,rl]=REAL[key];
    h+=`<tr><td>${key.replace(/\+/g,' + ')}</td>`+
       `<td class="n">${lk.a.toLocaleString()} km²<br>${lk.lv.toFixed(1)} m</td>`+
       `<td class="n">${ra.toLocaleString()} km²<br>${rl.toFixed(1)} m</td></tr>`;
  }
  document.getElementById('valid').innerHTML=h;
}

/* ---- transect section ---- */
const sc=document.getElementById('sec'), sx=sc.getContext('2d');
const T_SW=sc.width, T_SH=sc.height, T_PL=64, T_PR=18, T_PT=22, T_PB=40;
const Trans=RB_D.transect;
const DMAX=Trans.dist[Trans.dist.length-1];
function tx(d){return T_PL+d/DMAX*(T_SW-T_PL-T_PR);}
let T_YMIN=-260, T_YMAX=520;
function ty(z){return T_PT+(T_YMAX-z)/(T_YMAX-T_YMIN)*(T_SH-T_PT-T_PB);}
function drawSection(t){
  const fr=RB_F.find(f=>f.t===t)||RB_F[0];
  sx.fillStyle='#0f1317'; sx.fillRect(0,0,T_SW,T_SH);
  sx.strokeStyle='#232a31'; sx.lineWidth=1;
  sx.font='11px ui-monospace,Menlo,monospace'; sx.fillStyle='#71767c'; sx.textAlign='right';
  for(let z=-200;z<=500;z+=100){
    sx.beginPath(); sx.moveTo(T_PL,ty(z)); sx.lineTo(T_SW-T_PR,ty(z)); sx.stroke();
    sx.fillText(z+' m', T_PL-8, ty(z)+4);
  }
  sx.textAlign='center';
  for(let d=0;d<=DMAX;d+=100){ sx.fillText(Math.round(d)+' km', tx(d), T_SH-14); }
  sx.strokeStyle='#2f6f66'; sx.setLineDash([5,5]); sx.beginPath();
  sx.moveTo(T_PL,ty(fr.esl)); sx.lineTo(T_SW-T_PR,ty(fr.esl)); sx.stroke(); sx.setLineDash([]);
  sx.fillStyle='#2f6f66'; sx.textAlign='left';
  sx.fillText('sea level '+fr.esl.toFixed(0)+' m', T_PL+6, ty(fr.esl)-5);
  sx.strokeStyle='rgba(190,190,180,.35)'; sx.lineWidth=1.5; sx.beginPath();
  Trans.modern.forEach((z,i)=>{ i?sx.lineTo(tx(Trans.dist[i]),ty(z)):sx.moveTo(tx(Trans.dist[i]),ty(z)); });
  sx.stroke();
  sx.beginPath(); sx.moveTo(tx(0), ty(T_YMIN));
  fr.ground.forEach((z,i)=> sx.lineTo(tx(Trans.dist[i]), ty(z)));
  sx.lineTo(tx(DMAX), ty(T_YMIN)); sx.closePath();
  sx.fillStyle='#5d5a4c'; sx.fill();
  sx.strokeStyle='#141a1e'; sx.lineWidth=1.2; sx.stroke();
  sx.fillStyle='rgba(45,105,155,.85)';
  let run=null;
  for(let i=0;i<=fr.water.length;i++){
    const w=fr.water[i];
    if(w!=null && run===null) run=i;
    if((w==null||i===fr.water.length) && run!==null){
      sx.beginPath(); sx.moveTo(tx(Trans.dist[run]), ty(fr.water[run]));
      for(let k=run;k<i;k++) sx.lineTo(tx(Trans.dist[k]), ty(fr.water[k]));
      for(let k=i-1;k>=run;k--) sx.lineTo(tx(Trans.dist[k]), ty(fr.ground[k]));
      sx.closePath(); sx.fill();
      run=null;
    }
  }
  sx.fillStyle='rgba(228,240,245,.92)';
  run=null;
  for(let i=0;i<=fr.ice.length;i++){
    const on=i<fr.ice.length && fr.ice[i];
    if(on && run===null) run=i;
    if((!on||i===fr.ice.length) && run!==null){
      sx.beginPath(); sx.moveTo(tx(Trans.dist[run]), ty(T_YMAX));
      for(let k=run;k<i;k++) sx.lineTo(tx(Trans.dist[k]), ty(Math.max(fr.ground[k], fr.water[k]??-999)));
      sx.lineTo(tx(Trans.dist[i-1]), ty(T_YMAX)); sx.closePath(); sx.fill();
      run=null;
    }
  }
  sx.textAlign='center'; sx.font='italic 12px Georgia,serif';
  for(const m of Trans.marks){
    let bi=0,bd=1e9;
    Trans.lon.forEach((lo,i)=>{const d=Math.hypot(lo-m.lon,(Trans.lat[i]-m.lat)*1.4); if(d<bd){bd=d;bi=i;}});
    const X=tx(Trans.dist[bi]);
    sx.strokeStyle='rgba(200,120,100,.4)'; sx.lineWidth=1;
    sx.beginPath(); sx.moveTo(X,T_PT); sx.lineTo(X,T_SH-T_PB); sx.stroke();
    sx.fillStyle='#b9928a'; sx.fillText(m.name, X, T_PT-7);
  }
  sx.textAlign='left'; sx.font='12px ui-sans-serif,system-ui,sans-serif'; sx.fillStyle='#8e8b84';
  sx.fillText('south-west', T_PL, T_SH-26); sx.textAlign='right';
  sx.fillText('north-east', T_SW-T_PR, T_SH-26);
}

/* ---- readout (rebound site table) ---- */
function probe(t,lo,la,model){ const g=gridFor(t,model); return t===0?0:(giaAt(g,lo,la)); }
function buildReadout(t){
  const ro=document.getElementById('readout');
  const nm = RB_MODEL==='ice7' ? 'ICE-7G_NA (VM7)' : 'ICE-6G_C (VM5a)';
  const sites=[['Chicago outlet',-87.80,41.72],['Port Huron',-82.42,43.00],
               ['Sault',-84.35,46.50],['North Bay sill',-79.45,46.32],['Kingston',-76.46,44.23]];
  let h=`<p style="margin:0 0 10px">Rebound readout at solved timestep ${t/1000} ka, ${nm}.</p>`;
  h+='<table><tr><th>Site</th><th style="text-align:right">ICE-7G</th><th style="text-align:right">ICE-6G</th><th style="text-align:right">difference</th></tr>';
  for(const [n,lo,la] of sites){
    const a=probe(t,lo,la,'ice7'), b=probe(t,lo,la,'ice6');
    h+=`<tr><td>${n}</td><td class="n">${a==null?'—':a.toFixed(1)+' m'}</td>`+
       `<td class="n">${b==null?'—':b.toFixed(1)+' m'}</td>`+
       `<td class="n">${(a==null||b==null)?'—':((a-b>=0?'+':'')+(a-b).toFixed(1))+' m'}</td></tr>`;
  }
  const nb=probe(t,-79.45,46.32,RB_MODEL), ph=probe(t,-82.42,43.00,RB_MODEL);
  h+='</table>';
  if(nb!=null&&ph!=null)
    h+=`<p style="margin:10px 0 0">Differential between the North Bay and Port Huron sills at this date: <b>${(nb-ph).toFixed(0)} m</b>. `+
       `North Bay stands <b>${(nb-ph)>0?'lower':'higher'}</b> than Port Huron, which is what decides whether the upper lakes drain east down the Ottawa or south through the St. Clair.</p>`;
  ro.innerHTML=h;
}

/* ---- state ---- */
const yrEl=document.getElementById('yr'), snapEl=document.getElementById('snap'),
      stEl=document.getElementById('state'), confEl=document.getElementById('conf'),
      slider=document.getElementById('slider');

function nearest(y){
  let best=LK_F[0], bd=1e9;
  for(const fr of LK_F){ const d=Math.abs(fr.t-y); if(d<bd){bd=d;best=fr;} }
  return best;
}
function nearestT(y){
  let b=LK_F[0].t,bd=1e9;
  for(const f of LK_F){const d=Math.abs(f.t-y);if(d<bd){bd=d;b=f.t;}}
  return b;
}

function render(y){
  const {t,lkfr}=drawMap();
  yrEl.innerHTML = (y<50 ? 'present day' : Math.round(y).toLocaleString()+' <small>years before present</small>');
  snapEl.textContent = `solved timestep ${lkfr.t===0?'0':(lkfr.t/1000).toFixed(1)+' ka'} · `+
    `sea level ${lkfr.esl.toFixed(1)} m · ice covers ${lkfr.ice.length?'part of the domain':'nothing in the domain'}`;
  const rows=lkfr.lakes.filter(l=>l.a>=3000).sort((a,b)=>b.a-a.a).slice(0,7);
  let h='<tr><th>Water body</th><th style="text-align:right">level</th><th style="text-align:right">area</th><th style="text-align:right">sill</th></tr>';
  for(const lk of rows){
    h+=`<tr><td>${lk.b.length?lk.b.join(' + '):'unnamed basin'}</td>`+
       `<td class="n">${lk.lv.toFixed(1)} m</td>`+
       `<td class="n">${lk.a.toLocaleString()} km²</td>`+
       `<td class="n">${lk.sp?lk.sp[0].toFixed(2)+', '+lk.sp[1].toFixed(2):'—'}</td></tr>`;
  }
  stEl.innerHTML=h;
  confEl.textContent = 'Sill coordinates are found by locating where each lake touches ground that drains off the domain, then taking the lowest such cell. Lakes below the top six by area are not given a sill.';
  const X=hx(Math.min(14500,Math.max(0,y)));
  hyHead.setAttribute('x1',X); hyHead.setAttribute('x2',X); hyHd.setAttribute('cx',X);

  drawSection(t);
  buildReadout(t);
}

let cur=14500;
function set(y,fromSlider){
  cur=Math.max(0,Math.min(14500,y));
  if(!fromSlider) slider.value=Math.round((14500-cur)/100)*100;
  render(cur);
}
slider.addEventListener('input',()=>set(14500-(+slider.value),true));
document.getElementById('back').addEventListener('click',()=>{stop();set(cur+100);});
document.getElementById('fwd').addEventListener('click',()=>{stop();set(cur-100);});
let raf=null,last=0; const pb=document.getElementById('play');
function step(ts){ if(!last) last=ts; const dt=Math.min(64,ts-last); last=ts;
  let y=cur-dt*1.5; if(y<=0){set(0);stop();return;} set(y); raf=requestAnimationFrame(step); }
function play(){ if(cur<=0) set(14500); pb.textContent='Pause'; last=0; raf=requestAnimationFrame(step); }
function stop(){ if(raf) cancelAnimationFrame(raf); raf=null; pb.textContent='Play'; }
pb.addEventListener('click',()=> raf?stop():play());

function setMargin(m){
  LK_MARGIN=m; LK_D=LK_DS[m];
  LK_F=LK_D.frames.slice().sort((a,b)=>b.t-a.t);
  MODERN=LK_F[LK_F.length-1];
  pcache.clear(); buildHydro(); buildValid(); render(cur);
  document.getElementById('mgOpt').setAttribute('aria-pressed',String(m==='optimal'));
  document.getElementById('mgMax').setAttribute('aria-pressed',String(m==='max'));
}
document.getElementById('mgOpt').addEventListener('click',()=>setMargin('optimal'));
document.getElementById('mgMax').addEventListener('click',()=>setMargin('max'));

function setModel(m){
  RB_MODEL=m;
  document.getElementById('m6').setAttribute('aria-pressed',String(m==='ice6'));
  document.getElementById('m7').setAttribute('aria-pressed',String(m==='ice7'));
  render(cur);
}
document.getElementById('m6').addEventListener('click',()=>setModel('ice6'));
document.getElementById('m7').addEventListener('click',()=>setModel('ice7'));

for(const radio of document.querySelectorAll('input[name=base]')){
  radio.addEventListener('change',()=>{ if(radio.checked){ RB_MODE=radio.value; render(cur); } });
}
for(const id of ['lyLakes','lyIce','lySea','lyModern','lyAlg','lySpill','lyLabels','lyTransect']){
  document.getElementById(id).addEventListener('change',()=>render(cur));
}

buildHydro();
buildValid();
set(14500);
}

})();
