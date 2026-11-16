// Quantum Clock UI enhancements
export function installClockEnhancements() {
  // Collapsible advanced panels
  document.querySelectorAll('.box').forEach((box, i) => {
    if (i > 2) {
      const h = box.querySelector('h3');
      if (h) {
        h.style.cursor = 'pointer';
        h.onclick = () => box.querySelectorAll(':scope > *:not(h3)').forEach(el => el.style.display = el.style.display === 'none' ? '' : 'none');
      }
    }
  });

  // Matrix T-gate callout
  const callout = document.createElement('div');
  callout.className = 'box';
  callout.style.marginTop = '1rem';
  callout.innerHTML = '<h3>Matrix bridge</h3><p class="note">One forward tick of this clock (ωτ = π/4) is exactly the T gate of The Matrix. Export a receipt and open it in the Matrix lab to verify the handoff.</p>';
  const bridge = document.getElementById('bridge') || document.querySelector('[id*="matrix"]');
  if (bridge) bridge.prepend(callout);

  // Export + golden comparison button
  const btn = document.createElement('button');
  btn.className = 'btn';
  btn.textContent = 'Export receipt + golden Δ';
  btn.onclick = async () => {
    const s = window.clock?.state?.() || {};
    const payload = { schema: 'city.quantum.state-receipt.v1', tau: s.tau, phase: s.phaseTime, goldenDelta: Math.abs((s.phaseTime || 0) - (s.tau || 0)) };
    const blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'quantum-clock-receipt-golden.json';
    a.click();
  };
  const controls = document.querySelector('.row') || document.body;
  controls.appendChild(btn);

  // Contraction sparkline
  const spark = document.createElement('canvas');
  spark.width = 200; spark.height = 40;
  spark.style.marginTop = '8px';
  const ctx = spark.getContext('2d');
  function drawSpark() {
    ctx.clearRect(0, 0, 200, 40);
    ctx.strokeStyle = '#e3bf66';
    ctx.beginPath();
    for (let n = 0; n < 40; n++) {
      const y = 36 - 30 * Math.pow(15/17, n);
      if (n === 0) ctx.moveTo(n * 5, y); else ctx.lineTo(n * 5, y);
    }
    ctx.stroke();
  }
  drawSpark();
  const fig = document.getElementById('fig-contr') || document.body;
  fig.parentElement?.appendChild(spark);

  // Footer
  const footer = document.createElement('footer');
  footer.style.cssText = 'margin-top:2rem;padding:1rem;border-top:1px solid #273027;font:0.76rem IBM Plex Mono;color:#7f8375;display:flex;justify-content:space-between';
  footer.innerHTML = '<span>Numerical receipt · 48/48 precision · recomputed in browser</span><a href="https://nicholaskouns-create.github.io/E47-Kartekeya/cite/" style="color:#e3bf66" target="_blank">Open provenance</a>';
  document.body.appendChild(footer);
}
