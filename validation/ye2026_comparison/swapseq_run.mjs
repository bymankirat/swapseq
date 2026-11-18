// SwapSeq on one Ye et al. dataset: species relative abundance (%) x 1000 -> pseudo-counts; label = patient
import fs from 'fs';
import {analyse,annotate,distances,measurePower} from '../../method/core.mjs';
const [,,prof,meta,out,power]=process.argv;
const L=fs.readFileSync(prof,'utf8').split(/\r?\n/).filter(Boolean),h=L[0].split('\t').slice(1);
const M=new Map(fs.readFileSync(meta,'utf8').split(/\r?\n/).filter(Boolean).slice(1).map(l=>{const c=l.split('\t');return [c[0],c[1]]}));
const keep=h.map((id,j)=>M.has(id)?j:-1).filter(j=>j>=0);
const per=keep.map(()=>[]);const feats=[];
L.slice(1).forEach((l,fi)=>{const c=l.split('\t');feats.push(c[0]);keep.forEach((j,k)=>{const v=+c[j+1];if(v>0)per[k].push([fi,Math.max(1,Math.round(v*1000))])})});
const samples=keep.map((j,k)=>({id:h[j],idx:Int32Array.from(per[k].map(x=>x[0])),val:Float64Array.from(per[k].map(x=>x[1]))})).filter(s=>s.idx.length>0);
const lab=samples.map(s=>M.get(s.id));const data={feats,samples};
const Md=distances(data);const res=annotate(data,lab,analyse(data,lab,{},Md));
const pw=power==='1'?measurePower(data,lab,res,Md):null;
fs.writeFileSync(out,JSON.stringify({n:samples.length,subjects:new Set(lab).size,power:pw,rows:res.rows.filter(r=>r.level!=='ok').map(r=>({id:r.id,level:r.level,kind:r.kind,prob:r.prob,looksLike:r.looksLike,partner:r.partner,subject:r.subject,rare:r.rare?{own:r.rare.withOwn,other:r.rare.withOther}:null}))}));
