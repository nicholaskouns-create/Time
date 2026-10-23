import test from 'node:test';import assert from 'node:assert/strict';
import {issue,check,keys} from '../security/e47-signed-tick.mjs';
test('real clock signed receipt, replay and tamper rejection',()=>{
 const {privateKey,publicKey}=keys();
 const args={repository:'owner/repo',commit:'abc',run:'42',action:'deploy'};
 const receipt=issue({...args,privateKey,coordinateDt:0.1});
 assert.ok(receipt.tau>0);const seen=new Set();
 assert.equal(check(receipt,{...args,publicKey,seen}),true);
 assert.throws(()=>check(receipt,{...args,publicKey,seen}),/replay/);
 assert.throws(()=>check({...receipt,tau:999},{...args,publicKey,seen:new Set()}),/signature invalid/);
 assert.throws(()=>check(receipt,{...args,commit:'other',publicKey,seen:new Set()}),/scope mismatch/);
});
