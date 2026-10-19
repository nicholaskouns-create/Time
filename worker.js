// Runs the precision certificate off the main thread.
import { runCertificate } from './engine/certificate.js';
self.onmessage = async (e) => {
  try {
    const golden = await fetch(new URL('./data/golden.json', import.meta.url)).then(r => r.json()).catch(() => null);
    const cert = await runCertificate({ golden, environment: e.data?.environment || {}, onProgress: (k, title) => self.postMessage({ type: 'progress', k, title }) });
    self.postMessage({ type: 'done', cert });
  } catch (err) { self.postMessage({ type: 'error', message: String(err && err.message || err) }); }
};
