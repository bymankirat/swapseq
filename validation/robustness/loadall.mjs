import * as A from '../../atlas/code/load2.mjs';
export function all(){const o={...A.VALID,...A.NEW,...A.BOREN};for(const k of A.ploegNiches())o['ploeg_'+k.replace(/\W+/g,'_')]=()=>A.ploeg(k);return o}
export function get(key){let ds=all()[key]();if(!/relative/i.test(ds.type||'')){const keep=[];ds.data.samples.forEach((s,i)=>{if(s.val.reduce((a,b)=>a+b,0)>=1000)keep.push(i)});ds={...ds,data:{feats:ds.data.feats,samples:keep.map(i=>ds.data.samples[i])},lab:keep.map(i=>ds.lab[i]),time:ds.time&&keep.map(i=>ds.time[i])}}return ds}
if(process.argv[1].endsWith('loadall.mjs'))for(const k of Object.keys(all())){try{const d=get(k);console.log(k.padEnd(28),d.data.samples.length,new Set(d.lab).size)}catch(e){console.log(k.padEnd(28),'FAIL',e.message.slice(0,90))}}
