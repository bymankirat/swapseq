import {get as g2,all} from './loadall.mjs';import {MORE} from '../../atlas/code/load3.mjs';import {load as ye} from '../ye2026_comparison/loadye.mjs';
export const keysMore=Object.keys(MORE);
export function get(key){if(key.startsWith('ye:')){const {data,lab}=ye(key.slice(3));return {data,lab,type:'relative'}}
  if(MORE[key]){let ds=MORE[key]();const keep=[];ds.data.samples.forEach((s,i)=>{if(s.val.reduce((a,b)=>a+b,0)>=1000)keep.push(i)});return {...ds,data:{feats:ds.data.feats,samples:keep.map(i=>ds.data.samples[i])},lab:keep.map(i=>ds.lab[i]),house:ds.house&&keep.map(i=>ds.house[i]),time:keep.map(()=>0)}}
  return g2(key)}
