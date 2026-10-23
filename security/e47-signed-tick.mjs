// Real E47 NodeClock state, cryptographically signed and replay-checked.
// Requires Node 22. Each repo independently verifies with its configured public key.
import {NodeClock} from '../engine/sim.js';
import {generateKeyPairSync,sign,verify,createHash,randomUUID} from 'node:crypto';
import {readFileSync,writeFileSync,existsSync} from 'node:fs';
const stable=x=>JSON.stringify(x,Object.keys(x).sort());
const sha=x=>createHash('sha256').update(x).digest('hex');
export function issue({repository,commit,run,action,privateKey,clock=new NodeClock(),coordinateDt=1,nonce=randomUUID(),previous='0'.repeat(64)}){
 if(!repository||!commit||!run||!action||!privateKey||!(coordinateDt>0))throw Error('missing binding');
 clock.advance(coordinateDt);
 const state=clock.state();
 const body={schema:'e47-signed-tick-v1',repository,commit,run,action,nonce,previous,
  clock_id:'E47-NodeClock',tau:state.tau,coordinate_t:state.t,phase_time:state.phaseTime,
  lapse:state.lapse,observed_at_utc:new Date().toISOString()};
 const digest=sha(stable(body));
 return {...body,digest,signature:sign(null,Buffer.from(digest,'hex'),privateKey).toString('base64')};
}
export function check(receipt,{publicKey,repository,commit,run,action,seen=new Set(),maxAgeMs=300000,now=Date.now()}){
 if(!receipt||receipt.schema!=='e47-signed-tick-v1'||receipt.repository!==repository||receipt.commit!==commit||receipt.run!==run||receipt.action!==action)throw Error('scope mismatch');
 const {digest,signature,...body}=receipt;
 if(!Number.isFinite(body.tau)||body.tau<0||!Number.isFinite(body.phase_time)||!Number.isFinite(body.lapse))throw Error('invalid clock');
 const time=Date.parse(body.observed_at_utc);
 if(!Number.isFinite(time)||Math.abs(now-time)>maxAgeMs)throw Error('expired clock receipt');
 if(seen.has(body.nonce))throw Error('replay');
 if(sha(stable(body))!==digest||!verify(null,Buffer.from(digest,'hex'),publicKey,Buffer.from(signature,'base64')))throw Error('signature invalid');
 seen.add(body.nonce);return true;
}
export function keys(){return generateKeyPairSync('ed25519');}
