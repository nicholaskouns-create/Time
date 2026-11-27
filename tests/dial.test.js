// node --test tests/dial.test.js   (Node 22+, no dependencies)
// A Dial of Forty-Seven Hours: every check of the dial certificate, on the clock's own engine.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { runDialCertificate, DIAL_CERTIFICATE_ID } from '../dial/dial-certificate.js';

test(`${DIAL_CERTIFICATE_ID}: every check passes`, () => {
  const cert = runDialCertificate();
  for (const c of cert.checks) assert.ok(c.pass, `${c.id} ${c.title}: ${JSON.stringify(c.value)} residual ${c.residual}`);
  assert.equal(cert.passed, cert.total);
  assert.equal(cert.total, 12);
});
