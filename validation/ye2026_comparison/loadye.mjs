import fs from 'fs';
const D=(process.env.YE_DATA||'./Find-mislabels_Profiles_and_metadata')+'/';
export function load(name){const L=fs.readFileSync(D+name+'_profiles.txt','utf8').split(/\r?\n/).filter(Boolean),h=L[0].split('\t').slice(1);
  const M=new Map(fs.readFileSync(D+name+'_meta.txt','utf8').split(/\r?\n/).filter(Boolean).slice(1).map(l=>{const c=l.split('\t');return [c[0],c[1]]}));
  const keep=h.map((id,j)=>M.has(id)?j:-1).filter(j=>j>=0),per=keep.map(()=>[]),feats=[];
  L.slice(1).forEach((l,fi)=>{const c=l.split('\t');feats.push(c[0]);keep.forEach((j,k)=>{const v=+c[j+1];if(v>0)per[k].push([fi,Math.max(1,Math.round(v*1000))])})});
  const samples=keep.map((j,k)=>({id:h[j],idx:Int32Array.from(per[k].map(x=>x[0])),val:Float64Array.from(per[k].map(x=>x[1]))})).filter(s=>s.idx.length>0);
  return {data:{feats,samples},lab:samples.map(s=>M.get(s.id))}}
