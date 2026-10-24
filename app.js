const $ = (s) => document.querySelector(s);
const inserts = $('#inserts');
let lastReport = null;
let lastReports = [];
let loadedFastq = null;
let analysisBusy = false;

function cleanSequence(value) {
  return (value || '').toUpperCase().replace(/^>[^\n]*\n?/m, '').replace(/[^ACGTN]/g, '');
}
function parseSequence(value) {
  const raw = (value || '').trim();
  const header = raw.match(/^>([^\n\r]+)/)?.[1]?.trim() || '';
  return { name: header, sequence: cleanSequence(raw) };
}
function reverseComplement(seq) {
  return seq.split('').reverse().map(c => ({ A:'T', T:'A', C:'G', G:'C', N:'N' }[c] || 'N')).join('');
}
function parseFastq(value) {
  const lines = String(value || '').trim().split(/\r?\n/), reads = [];
  for (let i = 0; i + 3 < lines.length; i += 4) {
    if (!lines[i].startsWith('@') || lines[i + 2] !== '+') continue;
    const sequence = cleanSequence(lines[i + 1]);
    if (sequence) reads.push(sequence);
  }
  return reads;
}
function fastqConsensus(expected, reads) {
  const votes = Array.from({ length: expected.length }, () => ({ A:0, C:0, G:0, T:0, N:0 }));
  const k = 15, index = new Map();
  for (let i = 0; i <= expected.length - k; i++) { const key = expected.slice(i, i + k); if (!index.has(key)) index.set(key, []); index.get(key).push(i); }
  let mapped = 0;
  for (const original of reads) {
    const candidates = [original, reverseComplement(original)];
    let best = null;
    for (const read of candidates) {
      if (read.length < 12) continue;
      const anchors = [read.slice(0, k), read.slice(Math.max(0, Math.floor(read.length / 2) - Math.floor(k / 2)), Math.max(0, Math.floor(read.length / 2) - Math.floor(k / 2)) + k)];
      const candidatePositions = [...new Set(anchors.flatMap(anchor => index.get(anchor) || []).flatMap(pos => { const out=[]; for (let p=Math.max(0,pos-8); p<=Math.min(expected.length-read.length,pos+8); p++) out.push(p); return out; }))];
      for (const pos of candidatePositions) {
        let same = 0;
        for (let i = 0; i < read.length; i++) if (read[i] === expected[pos + i] || read[i] === 'N') same++;
        const identity = same / read.length;
        if (!best || identity > best.identity) best = { read, pos, identity };
      }
    }
    if (!best || best.identity < 0.75) continue;
    mapped++;
    for (let i = 0; i < best.read.length; i++) if (best.pos + i < votes.length) votes[best.pos + i][best.read[i]]++;
  }
  const consensus = votes.map(v => { const entries = Object.entries(v).sort((a,b) => b[1] - a[1]); return entries[0][1] ? entries[0][0] : 'N'; }).join('');
  const covered = votes.filter(v => Object.values(v).some(x => x > 0)).length;
  return { consensus, mapped, covered, total: reads.length };
}
function normalizeCircular(expected, observed) {
  if (!$('#circular').checked || expected.length !== observed.length) return observed;
  const anchor = expected.slice(0, Math.min(24, expected.length));
  const doubled = observed + observed;
  const rotation = doubled.indexOf(anchor);
  if (rotation > 0 && rotation < observed.length) return doubled.slice(rotation, rotation + observed.length);
  return observed;
}
function addInsert(name = `Insert ${inserts.children.length + 1}`, sequence = '') {
  const row = document.createElement('div'); row.className = 'insert-row';
  row.innerHTML = `<div class="insert-head"><label>Insert</label><div class="insert-actions"><button class="reverse" type="button">Reverse complement</button><button class="remove" type="button">Remove</button></div></div><div class="insert-fields"><input class="insert-name" value="${escapeHtml(name)}" placeholder="Name" /><textarea class="insert-seq" placeholder=">insert\nACGT...">${escapeHtml(sequence)}</textarea></div>`;
  row.querySelector('.remove').onclick = () => { row.remove(); renumber(); };
  row.querySelector('.reverse').onclick = () => {
    const field = row.querySelector('.insert-seq');
    const parsed = parseSequence(field.value);
    if (!parsed.sequence) return;
    field.value = `${parsed.name ? `>${parsed.name}\n` : ''}${reverseComplement(parsed.sequence)}`;
  };
  inserts.appendChild(row);
}
function renumber() { [...inserts.children].forEach((row, i) => { const input = row.querySelector('.insert-name'); if (!input.value) input.value = `Insert ${i + 1}`; }); }
function addClone(name = `Clone ${document.querySelectorAll('.clone-sample').length + 1}`, sequence = '') {
  const row = document.createElement('div'); row.className = 'clone-sample'; row._reads = null;
  row.innerHTML = `<div class="clone-sample-head"><input class="clone-name" value="${escapeHtml(name)}" placeholder="Sample name" /><button class="remove remove-clone" type="button">Remove</button></div><div class="file-row"><input class="clone-file" type="file" accept=".fa,.fasta,.fastq,.fq,.txt" /><button class="ghost clone-file-button" type="button">Choose FASTA / FASTQ</button></div><textarea class="clone-input" placeholder=">clone-01\nPaste FASTA/raw sequence, or upload FASTQ above">${escapeHtml(sequence)}</textarea><div class="sequence-meta"><span class="clone-length">${parseSequence(sequence).sequence.length.toLocaleString()} bp</span><span>Allowed: A C G T N</span></div>`;
  const textarea = row.querySelector('.clone-input'), length = row.querySelector('.clone-length'), file = row.querySelector('.clone-file');
  textarea.addEventListener('input', () => { row._reads = null; length.textContent = `${parseSequence(textarea.value).sequence.length.toLocaleString()} bp`; });
  row.querySelector('.clone-file-button').onclick = () => file.click();
  file.addEventListener('change', async () => {
    const selected = file.files[0]; if (!selected) return;
    const text = await selected.text();
    const isFastq = /\.(fastq|fq)$/i.test(selected.name) || text.trimStart().startsWith('@') && text.split(/\r?\n/)[2] === '+';
    if (isFastq) { row._reads = parseFastq(text); textarea.value = `FASTQ loaded: ${row._reads.length} reads\nClick Check cloning result to build a consensus.`; length.textContent = `${row._reads.length.toLocaleString()} reads`; }
    else { row._reads = null; textarea.value = text; length.textContent = `${parseSequence(text).sequence.length.toLocaleString()} bp`; }
  });
  row.querySelector('.remove-clone').onclick = () => { if (document.querySelectorAll('.clone-sample').length > 1) { row.remove(); } };
  $('#clones').appendChild(row);
  return row;
}
function escapeHtml(s) { return String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'}[c])); }
function readConstruct(observedForPlacement = '') {
  const rawBackbone = parseSequence($('#backboneInput').value).sequence;
  const removed = parseSequence($('#removedInput').value).sequence;
  const removedReverse = removed ? reverseComplement(removed) : '';
  const removedIndex = removed ? (rawBackbone.indexOf(removed) >= 0 ? rawBackbone.indexOf(removed) : rawBackbone.indexOf(removedReverse)) : -1;
  const strip = sequence => {
    if (!removed) return sequence;
    const forward = sequence.indexOf(removed);
    if (forward >= 0) return sequence.slice(0, forward) + sequence.slice(forward + removed.length);
    const reverse = sequence.indexOf(reverseComplement(removed));
    return reverse >= 0 ? sequence.slice(0, reverse) + sequence.slice(reverse + removed.length) : sequence;
  };
  const backbone = strip(rawBackbone);
  const parts = [...inserts.children].map((row, i) => ({ name: row.querySelector('.insert-name').value.trim() || `Insert ${i + 1}`, sequence: strip(parseSequence(row.querySelector('.insert-seq').value).sequence) })).filter(x => x.sequence);
  const direct = true;
  let inferredInsertAt = -1;
  if (!removed && !direct && observedForPlacement && parts.length) {
    for (const part of parts) {
      for (const target of [part.sequence, reverseComplement(part.sequence)]) {
        const observedIndex = observedForPlacement.indexOf(target);
        if (observedIndex < 0) continue;
        const leftFlank = observedForPlacement.slice(Math.max(0, observedIndex - 40), observedIndex);
        const leftIndex = leftFlank ? rawBackbone.lastIndexOf(leftFlank) : -1;
        if (leftIndex >= 0) { inferredInsertAt = leftIndex + leftFlank.length; break; }
        const rightFlank = observedForPlacement.slice(observedIndex + target.length, observedIndex + target.length + 40);
        const rightIndex = rightFlank ? rawBackbone.indexOf(rightFlank) : -1;
        if (rightIndex >= 0) { inferredInsertAt = rightIndex; break; }
      }
      if (inferredInsertAt >= 0) break;
    }
  }
  if (!removed && !direct && inferredInsertAt < 0 && observedForPlacement && parts.length) {
    const backboneAlignment = align(rawBackbone, observedForPlacement);
    const candidateBlocks = backboneAlignment ? findInsertedBlocks(backboneAlignment) : [];
    const matchedBlock = matchInsertedBlocks(candidateBlocks, parts).find(block => block.match && block.match.identity >= 0.6);
    if (matchedBlock) inferredInsertAt = matchedBlock.after;
  }
  const insertAt = removedIndex >= 0 ? removedIndex : inferredInsertAt >= 0 ? inferredInsertAt : backbone.length;
  const insertedSequence = parts.map(x => x.sequence).join('');
  const expected = direct ? backbone : backbone.slice(0, insertAt) + insertedSequence + backbone.slice(insertAt);
  let cursor = insertAt;
  const partRanges = parts.map(part => { const range = { start: cursor, end: cursor + part.sequence.length }; cursor = range.end; return range; });
  const backboneRanges = direct ? [{ start: 0, end: backbone.length }] : [{ start: 0, end: insertAt }, { start: insertAt + insertedSequence.length, end: expected.length }];
  return { backbone, parts, expected, direct, removed, insertAt, partRanges, backboneRanges, removalFound: !removed || removedIndex >= 0 };
}
function findBestMatch(expected, observed, start) {
  const window = expected.slice(start, start + Math.min(80, expected.length - start));
  if (!window) return { pos: -1, identity: 0 };
  let best = { pos: -1, identity: 0 };
  for (let p = 0; p <= observed.length - Math.min(12, window.length); p++) {
    let same = 0; const length = Math.min(window.length, observed.length - p);
    for (let i = 0; i < length; i++) if (window[i] === observed[p + i]) same++;
    const identity = length ? same / length : 0;
    if (identity > best.identity) best = { pos: p, identity };
  }
  return best;
}
function alignBanded(expected, observed) {
  // Banded Needleman-Wunsch: handles full plasmids around 10 kb while allowing
  // the small substitutions/indels expected from cloning and sequencing.
  const n = expected.length, m = observed.length, gap = -2, match = 2, mismatch = -3;
  const band = Math.min(1200, Math.max(250, Math.abs(n - m) + 120));
  if (Math.abs(n - m) > band) return null;
  const min = -2147483647, rows = [], lows = [];
  const scoreAt = (row, col) => {
    if (row < 0 || row >= rows.length || col < lows[row] || col > lows[row] + rows[row].length - 1) return min;
    return rows[row][col - lows[row]];
  };
  for (let i = 0; i <= n; i++) {
    const lo = Math.max(0, i - band), hi = Math.min(m, i + band);
    if (hi < lo) return null;
    const row = new Int32Array(hi - lo + 1); row.fill(min); lows.push(lo); rows.push(row);
    for (let j = lo; j <= hi; j++) {
      let value = min;
      if (i === 0 && j === 0) value = 0;
      else {
        if (i && j) value = Math.max(value, scoreAt(i - 1, j - 1) + (expected[i - 1] === observed[j - 1] ? match : mismatch));
        if (i) value = Math.max(value, scoreAt(i - 1, j) + gap);
        if (j) value = Math.max(value, scoreAt(i, j - 1) + gap);
      }
      row[j - lo] = value;
    }
  }
  if (scoreAt(n, m) === min) return null;
  let a='', b='', i=n, j=m;
  while (i || j) {
    const current = scoreAt(i, j);
    const diagonal = i && j ? scoreAt(i - 1, j - 1) + (expected[i - 1] === observed[j - 1] ? match : mismatch) : min;
    if (i && j && current === diagonal) { a=expected[i-1]+a; b=observed[j-1]+b; i--; j--; }
    else if (i && current === scoreAt(i - 1, j) + gap) { a=expected[i-1]+a; b='-'+b; i--; }
    else { a='-'+a; b=observed[j-1]+b; j--; }
  }
  return { expected: a, observed: b };
}
function alignByStructuralOffset(expected, observed) {
  const k = 25, index = new Map(), anchors = [];
  for (let i = 0; i <= expected.length - k; i++) { const key = expected.slice(i, i + k); if (!index.has(key)) index.set(key, []); index.get(key).push(i); }
  for (let j = 0; j <= observed.length - k; j++) { const hits = index.get(observed.slice(j, j + k)); if (hits && hits.length === 1) anchors.push({ e: hits[0], o: j, offset: j - hits[0] }); }
  anchors.sort((a, b) => a.e - b.e);
  for (let i = 1; i < anchors.length; i++) {
    const before = anchors[i - 1], after = anchors[i], jump = after.offset - before.offset;
    if (jump < 30 || jump > 2000 || after.e - before.e < k || after.o - before.o < k) continue;
    const beforeEnd = before.e + k, afterStart = after.e, beforeObservedEnd = before.o + k, afterObservedStart = after.o;
    const left = alignBanded(expected.slice(0, beforeEnd), observed.slice(0, beforeObservedEnd));
    const right = alignBanded(expected.slice(afterStart), observed.slice(afterObservedStart));
    if (!left || !right) continue;
    const expectedGap = expected.slice(beforeEnd, afterStart), observedGap = observed.slice(beforeObservedEnd, afterObservedStart);
    if (observedGap.length - expectedGap.length !== jump) continue;
    const extra = observedGap.length - expectedGap.length;
    const middleExpected = '-'.repeat(extra) + expectedGap;
    const middleObserved = observedGap;
    return { expected: left.expected + middleExpected + right.expected, observed: left.observed + middleObserved + right.observed };
  }
  return null;
}
function align(expected, observed) {
  if (Math.max(expected.length, observed.length) < 2000) return alignBanded(expected, observed);
  const structural = alignByStructuralOffset(expected, observed);
  if (structural) return structural;
  const k=21, expectedIndex=new Map();
  for (let i=0; i<=expected.length-k; i++) { const key=expected.slice(i,i+k); if (!expectedIndex.has(key)) expectedIndex.set(key, []); expectedIndex.get(key).push(i); }
  const anchors=[];
  for (let j=0; j<=observed.length-k; j++) { const hits=expectedIndex.get(observed.slice(j,j+k)); if (hits && hits.length===1) anchors.push({e:hits[0],o:j}); }
  const selected=[]; let lastE=-1,lastO=-1;
  for (const anchor of anchors) {
    if (anchor.e<=lastE || anchor.o<=lastO) continue;
    if (!selected.length || anchor.e-lastE>=80 || anchor.o-lastO>=80) { selected.push(anchor); lastE=anchor.e; lastO=anchor.o; }
  }
  if (selected.length < 4) return alignBanded(expected, observed);
  let a='', b='', prevE=0, prevO=0;
  for (const anchor of selected) {
    const left=alignBanded(expected.slice(prevE,anchor.e), observed.slice(prevO,anchor.o));
    if (!left) return alignBanded(expected, observed);
    a+=left.expected; b+=left.observed; a+=expected.slice(anchor.e,anchor.e+k); b+=observed.slice(anchor.o,anchor.o+k); prevE=anchor.e+k; prevO=anchor.o+k;
  }
  const tail=alignBanded(expected.slice(prevE), observed.slice(prevO));
  if (!tail) return alignBanded(expected, observed);
  return { expected:a+tail.expected, observed:b+tail.observed };
}
function analyze(event) {
  event?.preventDefault();
  if (analysisBusy) return;
  analysisBusy = true;
  const button = $('#analyze');
  button.disabled = true;
  button.classList.add('loading');
  button.innerHTML = 'Checking…';
  try {
    const rows = [...document.querySelectorAll('.clone-sample')];
    if (!rows.length) return alert('Add at least one cloned plasmid sample.');
    const reports = [];
    for (const row of rows) {
      const sampleName = row.querySelector('.clone-name').value.trim() || `Clone ${reports.length + 1}`;
      let observed = parseSequence(row.querySelector('.clone-input').value).sequence, fastqStats = null;
      const construct = readConstruct(observed);
      if (row._reads) { fastqStats = fastqConsensus(construct.expected, row._reads); observed = fastqStats.consensus; }
      observed = normalizeCircular(construct.expected, observed);
      if (!construct.backbone || !construct.expected || !observed) return alert(`Please provide a valid sequence for ${sampleName}.`);
      const aln = align(construct.expected, observed);
      if (!aln) return alert(`${sampleName} differs by more than the supported alignment band.`);
      const diffs = []; let expPos = 0;
      for (let i = 0; i < aln.expected.length; i++) { const e=aln.expected[i], o=aln.observed[i]; if (e !== o) diffs.push({ position: expPos + 1, expected: e, observed: o, type: e==='-'?'Insertion':o==='-'?'Deletion':'Substitution' }); if (e !== '-') expPos++; }
      const regions = [regionResult('Backbone', construct.backbone, aln, construct.backboneRanges)];
      if (!construct.direct) construct.parts.forEach((part, i) => regions.push(regionResult(part.name, part.sequence, aln, [construct.partRanges[i]])));
      const insertions = matchInsertedBlocks(findInsertedBlocks(aln), construct.parts);
      reports.push({ sampleName, construct, observed, diffs, regions, insertions, alignment: aln, matches: diffs.filter(d => d.type === 'Substitution').length, fastqStats });
    }
    renderMultiReport(reports);
  } catch (error) {
    console.error(error);
    alert(`Analysis could not be completed: ${error.message || error}`);
  } finally {
    analysisBusy = false;
    button.disabled = false;
    button.classList.remove('loading');
    button.innerHTML = 'Check again <span>→</span>';
  }
}
function regionResult(name, seq, aln, ranges) {
  let expectedPos = 0, expectedChars = 0, observedChars = 0, errors = 0;
  for (let i = 0; i < aln.expected.length; i++) {
    const e = aln.expected[i], o = aln.observed[i];
    if (e !== '-') {
      const inRegion = ranges.some(range => expectedPos >= range.start && expectedPos < range.end);
      if (inRegion) { expectedChars++; if (o !== '-') { observedChars++; if (e !== o) errors++; } else errors++; }
      expectedPos++;
    }
  }
  const exact = expectedChars === seq.length && observedChars === seq.length && errors === 0;
  const covered = observedChars >= Math.max(1, expectedChars * .95);
  return { name, expected: expectedChars, observed: observedChars, status: exact ? 'pass' : covered ? 'warn' : 'fail' };
}
function compareSequences(a, b) {
  const aln = align(a, b);
  if (!aln) return { identity: 0, differences: Infinity };
  let matches=0, differences=0, compared=0;
  for (let i=0; i<aln.expected.length; i++) { if (aln.expected[i] !== '-' || aln.observed[i] !== '-') compared++; if (aln.expected[i] !== aln.observed[i]) differences++; else matches++; }
  return { identity: compared ? matches / compared : 0, differences };
}
function localCompareSequences(a, b) {
  const n = a.length, m = b.length;
  if (!n || !m) return { identity: 0, alignedLength: 0, predictedCoverage: 0, insertCoverage: 0, score: 0 };
  const rowWidth = m + 1, scores = new Int32Array(rowWidth), traces = new Uint8Array((n + 1) * rowWidth);
  let bestScore = 0, bestI = 0, bestJ = 0;
  for (let i = 1; i <= n; i++) {
    let diagonal = 0;
    for (let j = 1; j <= m; j++) {
      const up = scores[j] - 2, left = scores[j - 1] - 2, diag = diagonal + (a[i - 1] === b[j - 1] ? 2 : -3);
      diagonal = scores[j];
      let score = 0, trace = 0;
      if (diag > score) { score = diag; trace = 1; }
      if (up > score) { score = up; trace = 2; }
      if (left > score) { score = left; trace = 3; }
      scores[j] = score; traces[i * rowWidth + j] = trace;
      if (score > bestScore) { bestScore = score; bestI = i; bestJ = j; }
    }
  }
  let i = bestI, j = bestJ, matches = 0, alignedLength = 0, alignedA = '', alignedB = '';
  while (i > 0 && j > 0) {
    const trace = traces[i * rowWidth + j];
    if (!trace) break;
    if (trace === 1) { if (a[i - 1] === b[j - 1]) matches++; alignedLength++; alignedA=a[i-1]+alignedA; alignedB=b[j-1]+alignedB; i--; j--; }
    else if (trace === 2) { alignedLength++; alignedA=a[i-1]+alignedA; alignedB='-'+alignedB; i--; }
    else { alignedLength++; alignedA='-'+alignedA; alignedB=b[j-1]+alignedB; j--; }
  }
  return { identity: alignedLength ? matches / alignedLength : 0, alignedLength, predictedCoverage: alignedLength / n, insertCoverage: alignedLength / m, score: bestScore, alignedA, alignedB, predictedStart: i, insertStart: j };
}
function matchInsertedBlocks(blocks, parts) {
  return blocks.map(block => {
    if (!parts.length) return { ...block, match:null };
    let best=null;
    for (const part of parts) for (const orientation of ['forward','reverse-complement']) {
      const target=orientation==='forward' ? part.sequence : reverseComplement(part.sequence), score=localCompareSequences(block.insertedSequence || block.sequence,target);
      if (!best || score.score > best.score) best={ name:part.name, orientation, targetSequence:target, ...score };
    }
    return { ...block, match:best };
  });
}
function renderMultiReport(reports) {
  lastReports = reports;
  const tabs = $('#sampleTabs');
  tabs.classList.toggle('hidden', reports.length < 2);
  tabs.innerHTML = reports.map((report, i) => `<button type="button" class="sample-tab ${sampleReportStatus(report)}${i === 0 ? ' active' : ''}" data-sample-index="${i}">${escapeHtml(report.sampleName)}</button>`).join('');
  tabs.querySelectorAll('.sample-tab').forEach(tab => tab.onclick = () => {
    tabs.querySelectorAll('.sample-tab').forEach(x => x.classList.remove('active'));
    tab.classList.add('active');
    renderReport(lastReports[Number(tab.dataset.sampleIndex)]);
  });
  lastReport = reports[0];
  renderReport(lastReport);
}
function sampleReportStatus(report) {
  const matched = report.insertions.some(x => x.match && x.match.identity >= .95 && x.match.predictedCoverage >= .85 && x.match.insertCoverage >= .85);
  const partial = report.insertions.some(x => x.match && x.match.identity >= .8);
  if (report.construct.direct) return matched ? 'match' : partial ? 'partial' : 'none';
  if (report.diffs.length === 0 && report.observed.length === report.construct.expected.length) return 'match';
  return partial ? 'partial' : 'none';
}
function renderReport(r) {
  $('#results').classList.remove('hidden');
  const pass = r.diffs.length === 0 && r.observed.length === r.construct.expected.length;
  const statusTitle = r.construct.direct ? (r.insertions.length ? 'Internal inserted sequence detected' : (pass ? 'Recovered plasmid matches the backbone' : 'Differences from backbone detected')) : (pass ? 'Cloning result matches the reference' : 'Cloning result needs review');
  $('#resultTitle').textContent = `${r.sampleName ? `${r.sampleName} · ` : ''}${statusTitle}`;
  $('#resultTitle').className = pass ? 'good' : 'bad';
  const coverage = r.fastqStats ? `${(r.fastqStats.covered / r.construct.expected.length * 100).toFixed(1)}%` : 'n/a';
  $('#summary').innerHTML = metric(pass ? 'PASS' : 'REVIEW', 'Overall status', pass ? 'good' : 'warn') + metric(`${r.observed.length.toLocaleString()} bp`, 'Consensus length') + metric(`${r.construct.expected.length.toLocaleString()} bp`, 'Expected length') + metric(r.fastqStats ? coverage : r.diffs.length, r.fastqStats ? 'Read coverage' : 'Differences', r.fastqStats ? (r.fastqStats.covered / r.construct.expected.length > .95 ? 'good' : 'warn') : (r.diffs.length ? 'bad' : 'good'));
  const regionRows = r.regions.map(x => `<div class="region"><div class="region-name">${escapeHtml(x.name)}</div><div class="region-detail">Expected ${x.expected.toLocaleString()} bp · observed ${x.observed.toLocaleString()} bp</div><span class="badge ${x.status}">${x.status === 'pass' ? 'PASS' : x.status === 'warn' ? 'REVIEW' : 'FAIL'}</span></div>`);
  const insertionRows = r.insertions.map((x, i) => { const m=x.match, matched=m && m.identity >= .95 && m.predictedCoverage >= .85 && m.insertCoverage >= .85, partial=m && m.identity >= .8, status=matched?'match':partial?'partial':'none', label=m?.name || `Inserted sequence ${i + 1}`; const detail=`${x.length.toLocaleString()} bp inserted after backbone position ${x.after.toLocaleString()}` + (m ? ` · ${escapeHtml(m.name)} ${(m.identity*100).toFixed(1)}% identity · ${Math.round(m.predictedCoverage*100)}% predicted / ${Math.round(m.insertCoverage*100)}% supplied · ${m.orientation}` : ' · no supplied insert to compare'); return `<div class="region"><div class="region-name insert-status-name ${status}">${escapeHtml(label)}</div><div class="region-detail">${detail}</div><span class="badge ${matched?'pass':partial?'warn':'fail'}">${matched?'MATCH':partial?'PARTIAL':'INSERT'}</span></div>`; });
  $('#regionsPanel').innerHTML = regionRows.concat(insertionRows).join('') || '<p class="mut">No regions available.</p>';
  $('#differencesPanel').innerHTML = r.diffs.length ? r.diffs.slice(0,250).map(d => `<div class="diff"><span><strong>${d.type}</strong> at expected position ${d.position.toLocaleString()}</span><code>${d.expected} → ${d.observed}</code></div>`).join('') + (r.diffs.length > 250 ? '<p class="mut">Showing first 250 differences.</p>' : '') : '<p class="good">No sequence differences detected.</p>';
  $('#sequencePanel').innerHTML = `${detailedAlignmentViewer(r)}<p class="mut">Full alignment text. Mismatches are marked with ^.</p><div class="sequence-box">EXPECTED  ${escapeHtml(r.alignment.expected)}\nOBSERVED  ${escapeHtml(r.alignment.observed)}</div>`;
  document.querySelectorAll('.tab').forEach(t => t.onclick = () => { document.querySelectorAll('.tab').forEach(x => x.classList.remove('active')); t.classList.add('active'); ['regions','differences','sequence'].forEach(n => $(`#${n}Panel`).classList.toggle('hidden', t.dataset.tab !== n)); });
  const copyButton = document.querySelector('.copy-sequence');
  if (copyButton) copyButton.onclick = async () => { await navigator.clipboard?.writeText(r.observed); copyButton.textContent='Copied'; setTimeout(() => copyButton.textContent='Copy sequence', 1200); };
  const alignmentScroller=document.querySelector('.sequence-scroll');
  const jumpToColumn=(column)=>{ if(!alignmentScroller) return; const ratio=Math.max(0,Math.min(1,Number(column)/Math.max(1,r.alignment.expected.length))); alignmentScroller.scrollLeft=ratio*alignmentScroller.scrollWidth-alignmentScroller.clientWidth/2; };
  document.querySelectorAll('.mismatch-marker,.difference-region').forEach(el=>el.addEventListener('click',()=>jumpToColumn(el.dataset.column)));
  $('#results').scrollIntoView({ behavior: 'smooth', block: 'start' });
}
function findInsertedBlocks(aln) {
  const blocks=[]; let active=null, expectedPos=0;
  for (let i=0; i<aln.expected.length; i++) {
    const e=aln.expected[i], o=aln.observed[i];
    if (e !== '-') expectedPos++;
    if (e === '-' && o !== '-') { if (!active) active={ after:expectedPos, sequence:'', insertedSequence:'', startColumn:i }; active.sequence += o; active.insertedSequence += o; }
    else if (active) { blocks.push({ after:active.after, length:active.insertedSequence.length, sequence:active.sequence, insertedSequence:active.insertedSequence, startColumn:active.startColumn, endColumn:i }); active=null; }
  }
  if (active) blocks.push({ after:active.after, length:active.insertedSequence.length, sequence:active.sequence, insertedSequence:active.insertedSequence, startColumn:active.startColumn, endColumn:aln.observed.length });
  const merged=[];
  for (const block of blocks) {
    const previous=merged[merged.length-1];
    if (previous && block.startColumn-previous.endColumn <= 100 && block.after-previous.after <= 100) {
      previous.endColumn=block.endColumn;
      previous.sequence=aln.observed.slice(previous.startColumn, previous.endColumn).replace(/-/g,'');
      previous.insertedSequence += block.insertedSequence;
      previous.length += block.length;
    } else merged.push(block);
  }
  return merged;
}
function alignmentDiagram(r) {
  const exp = r.alignment.expected, obs = r.alignment.observed, width = 1000, left = 48, right = 18, usable = width - left - right, step = usable / Math.max(1, exp.length);
  const mismatchPositions = [], columns = [];
  for (let i = 0; i < exp.length; i++) { if (exp[i] !== obs[i]) mismatchPositions.push(i); if (exp[i] !== '-') columns.push({ i, base: exp[i] }); }
  const segments = [];
  let cursor = 0;
  for (const part of r.construct.parts) { const start = r.construct.backbone.length + cursor, end = start + part.sequence.length; segments.push({ start, end, name: part.name }); cursor += part.sequence.length; }
  const segmentRects = segments.map((s, i) => { const x=left+s.start*step, w=Math.max(2,(s.end-s.start)*step); return `<rect x="${x.toFixed(2)}" y="38" width="${w.toFixed(2)}" height="20" rx="3" fill="${i%2?'#d9eee3':'#c4e4d2'}"/><text x="${(x+w/2).toFixed(2)}" y="33" text-anchor="middle" class="svg-label">${escapeHtml(s.name)}</text>`; }).join('');
  const mismatchRects = mismatchPositions.map(pos => { const x=left+pos*step; const e=exp[pos]==='-'?'insert':exp[pos], o=obs[pos]==='-'?'deleted':obs[pos]; return `<rect x="${Math.max(left,x-1).toFixed(2)}" y="70" width="${Math.max(2,step).toFixed(2)}" height="54" rx="2" fill="#e46b6b"><title>Alignment column ${pos+1}: ${e} → ${o}</title></rect>`; }).join('');
  const ticks = [0, .25, .5, .75, 1].map(p => { const x=left+usable*p; return `<line x1="${x}" y1="132" x2="${x}" y2="138" stroke="#9aa9a2"/><text x="${x}" y="153" text-anchor="middle" class="svg-axis">${Math.round(exp.length*p).toLocaleString()}</text>`; }).join('');
  return `<div class="alignment-figure"><p class="alignment-caption"><strong>Alignment map</strong> · insert regions in green · mismatches, insertions, or deletions in red</p><svg viewBox="0 0 ${width} 170" role="img" aria-label="Visual alignment map"><style>.svg-label{font:11px system-ui;fill:#35634f}.svg-axis{font:10px system-ui;fill:#6d7a74}</style><text x="4" y="53" class="svg-axis">REF</text><rect x="${left}" y="38" width="${usable}" height="20" rx="3" fill="#edf2ef"/>${segmentRects}<text x="4" y="86" class="svg-axis">OBS</text><rect x="${left}" y="70" width="${usable}" height="20" rx="3" fill="#eaf0ed"/>${mismatchRects}<line x1="${left}" y1="132" x2="${left+usable}" y2="132" stroke="#cbd6d0"/>${ticks}</svg><div class="alignment-legend"><span class="legend-item"><i class="legend-swatch" style="background:#c4e4d2"></i>Expected insert</span><span class="legend-item"><i class="legend-swatch" style="background:#e46b6b"></i>Mismatch / indel</span><span class="legend-item"><i class="legend-swatch" style="background:#eaf0ed"></i>Matching sequence</span></div></div>`;
}
function detailedAlignmentViewer(r) {
  const exp = r.alignment.expected, obs = r.alignment.observed, width = 1200, left = 62, right = 22, usable = width-left-right;
  const differences = [], expCoords = [];
  let coord = 0;
  for (let i=0; i<exp.length; i++) { if (exp[i] !== '-') coord++; expCoords[i]=coord; if (exp[i] !== obs[i]) differences.push({ column:i, position:exp[i] === '-' ? coord : expCoords[i], expected:exp[i], observed:obs[i] }); }
  const total = Math.max(1, r.construct.expected.length), scale = usable/total;
  const insertions = r.construct.direct ? [] : r.construct.parts.map((part, i) => { const range=r.construct.partRanges[i]; return {name:part.name,start:range.start,end:range.end,color:i%2?'#b9de23':'#8bb400'}; });
  const featureRects = insertions.map(x => { const px=left+x.start*scale, w=Math.max(9,(x.end-x.start)*scale); return `<g><rect x="${px.toFixed(1)}" y="45" width="${w.toFixed(1)}" height="26" rx="6" fill="${x.color}" stroke="#234100" stroke-width="2"/><text x="${(px+w/2).toFixed(1)}" y="39" text-anchor="middle" class="map-label">${escapeHtml(x.name)}</text></g>`; }).join('');
  const mismatchMarkers = differences.map((d,i) => { const p=Math.max(0,Math.min(total,d.position-1)), x=left+p*scale; return `<g class="mismatch-marker" data-column="${d.column}" tabindex="0"><line x1="${x.toFixed(1)}" y1="78" x2="${x.toFixed(1)}" y2="111" stroke="#e33b2f" stroke-width="${i===0?3:2}"/><circle cx="${x.toFixed(1)}" cy="116" r="5" fill="#e33b2f"><title>${d.expected==='-'?'Insertion':d.observed==='-'?'Deletion':'Mismatch'} at ${d.position}: ${d.expected==='-'?'—':d.expected} → ${d.observed==='-'?'—':d.observed}</title></circle></g>`; }).join('');
  const sampleRow = alignmentRow(obs, exp, 0, obs.length), refRow = alignmentRow(exp, obs, 0, exp.length);
  const differenceRegions = groupDifferenceRegions(exp, obs, expCoords);
  const ticks = [0,.25,.5,.75,1].map(p => { const x=left+usable*p; return `<line x1="${x}" y1="131" x2="${x}" y2="137" stroke="#aab6b0"/><text x="${x}" y="151" text-anchor="middle" class="map-axis">${Math.round(total*p).toLocaleString()}</text>`; }).join('');
  const mismatchText = differences.length ? `${differences.length} difference${differences.length===1?'':'s'} detected` : 'No differences detected';
  const regionButtons=differenceRegions.map((x,i)=>`<button type="button" class="difference-region" data-column="${x.startColumn}"><span>${x.kind}</span><strong>${x.start.toLocaleString()}–${x.end.toLocaleString()}</strong><small>${x.length.toLocaleString()} bp</small></button>`).join('');
  const bestInsert = r.insertions.map(block => ({ block, match: block.match })).filter(x => x.match).sort((a,b) => b.match.score - a.match.score)[0];
  const bestInsertRow = bestInsert ? (() => {
    const match = bestInsert.match;
    const fullOverlay = Array(exp.length).fill(' '), matchedOverlay = Array(exp.length).fill(false), guideOverlay = Array(exp.length).fill('<span class="insert-guide blank-guide"></span>');
    const floatingBases = [];
    let firstMappedColumn = exp.length, lastMappedColumn = -1;
    let observedIndex = match.predictedStart;
    for (let i = 0; i < match.alignedA.length; i++) {
      const observedBase = match.alignedA[i], insertBase = match.alignedB[i];
      if (observedBase === '-') { if (insertBase !== '-') floatingBases.push(insertBase); continue; }
      const column = bestInsert.block.startColumn + observedIndex;
      if (column >= 0 && column < exp.length) {
        if (insertBase !== '-') fullOverlay[column] = insertBase;
        else fullOverlay[column] = '·';
        firstMappedColumn = Math.min(firstMappedColumn, column);
        lastMappedColumn = Math.max(lastMappedColumn, column);
        if (insertBase !== '-' && observedBase === insertBase) { matchedOverlay[column] = true; guideOverlay[column] = '<span class="insert-guide"></span>'; }
      }
      observedIndex++;
    }
    const place = (bases, start, step) => {
      let position = start;
      for (const base of bases) {
        while (position >= 0 && position < exp.length && fullOverlay[position] !== ' ') position += step;
        if (position < 0 || position >= exp.length) break;
        fullOverlay[position] = base;
        position += step;
      }
    };
    const prefix = [...match.targetSequence.slice(0, match.insertStart)];
    const suffix = [...match.targetSequence.slice(match.insertStart + match.alignedB.replace(/-/g, '').length)];
    const anchor = firstMappedColumn <= lastMappedColumn ? firstMappedColumn : Math.max(0, bestInsert.block.startColumn + match.predictedStart);
    place(prefix, anchor - prefix.length, 1);
    place(floatingBases, anchor - floatingBases.length - prefix.length, 1);
    place(suffix, lastMappedColumn + 1, 1);
    const fullOverlayHtml = [...fullOverlay].map((base, i) => `<span class="insert-base${matchedOverlay[i] ? ' insert-match-base' : ''}">${escapeHtml(base)}</span>`).join('');
    return `<div class="sequence-line insert-alignment-line" title="${escapeHtml(match.name)} · ${Math.round(match.identity*100)}% local identity · full supplied sequence shown"><span class="sequence-label insert-label">Best insert</span><code>${fullOverlayHtml}</code></div><div class="sequence-line insert-guide-line"><span class="sequence-label"></span><code>${guideOverlay.join('')}</code></div>`;
  })() : '';
  return `<div class="reference-view"><div class="reference-view-head"><div><div class="reference-kicker">MISMATCH TO REFERENCE</div><h3>Recovered plasmid vs reference</h3></div><div class="reference-actions"><button class="copy-sequence" type="button">Copy sequence</button><strong>${r.observed.length.toLocaleString()} bp</strong></div></div><div class="reference-map"><div class="map-title"><span>Reference features</span><span class="${differences.length?'mismatch-count':'match-count'}">${mismatchText}</span></div><svg viewBox="0 0 ${width} 165" role="img" aria-label="Reference feature map with mismatch markers"><style>.map-label{font:12px system-ui;fill:#355522;font-weight:700}.map-axis{font:11px ui-monospace,monospace;fill:#69766f}.mismatch-marker{cursor:pointer}</style><line x1="${left}" y1="58" x2="${left+usable}" y2="58" stroke="#1e241f" stroke-width="3"/><circle cx="${left}" cy="58" r="7" fill="#fff" stroke="#1e241f" stroke-width="3"/>${featureRects}<rect x="${left}" y="84" width="${usable}" height="18" rx="4" fill="#edf0ed"/><text x="8" y="97" class="map-axis">SAMPLE</text>${mismatchMarkers}<line x1="${left}" y1="131" x2="${left+usable}" y2="131" stroke="#d2d9d4"/>${ticks}</svg><div class="map-legend"><span><i class="map-swatch insert-swatch"></i>insert/reference feature</span><span><i class="map-swatch mismatch-swatch"></i>mismatch or indel</span><span class="map-help">Click markers or regions to jump</span></div>${regionButtons?`<div class="difference-regions"><div class="difference-regions-title">Difference regions</div>${regionButtons}</div>`:''}</div><div class="sequence-window"><div class="window-head"><span>Full sequence alignment</span><span>${exp.length.toLocaleString()} alignment columns · scroll horizontally</span></div><div class="sequence-scroll">${bestInsertRow}<div class="sequence-line"><span class="sequence-label sample-label">Sample</span><code>${sampleRow}</code></div><div class="sequence-line"><span class="sequence-label reference-label">Reference</span><code>${refRow}</code></div></div><div class="window-legend"><span class="sample-key">Sample</span><span class="reference-key">Reference</span><span class="insert-key">⋮ = matching insert base</span><span>red = mismatch / insertion / deletion</span></div></div></div>`;
}
function groupDifferenceRegions(exp, obs, coords) {
  const regions=[]; let active=null; const mergeGap=100;
  for (let i=0; i<exp.length; i++) {
    if (exp[i]===obs[i]) continue;
    const kind=exp[i]==='-'?'Insertion':obs[i]==='-'?'Deletion':'Mismatch';
    if(!active || i-active.lastColumn>mergeGap){if(active)regions.push(active);active={kind,start:coords[i]||1,end:coords[i]||1,startColumn:i,lastColumn:i,length:1};}
    else {if(active.kind!==kind) active.kind=active.kind==='Insertion'||kind==='Insertion'?'Insertion':'Mixed'; active.end=coords[i]||active.end;active.lastColumn=i;active.length++;}
  }
  if(active)regions.push(active); return regions;
}
function alignmentRow(primary, comparison, start, end) {
  let html='';
  for (let i=start; i<end; i++) { const mismatch=primary[i]!==comparison[i], char=primary[i]==='-'?'·':primary[i]; html += `<span class="${mismatch?'sequence-mismatch':''}" title="Alignment column ${i+1}">${char}</span>`; }
  return html;
}
function fixedSequenceChars(sequence, className = '') {
  return [...sequence].map((char, i) => `<span class="${className}" title="Alignment column ${i + 1}">${escapeHtml(char)}</span>`).join('');
}
function metric(value, label, cls='') { return `<div class="metric"><div class="value ${cls}">${value}</div><div class="label">${label}</div></div>`; }
function loadFile(input, target) { input.addEventListener('change', async () => { const f=input.files[0]; if (!f) return; const text=await f.text(); const isFastq=/\.(fastq|fq)$/i.test(f.name) || text.trimStart().startsWith('@') && text.split(/\r?\n/)[2] === '+'; if (target.id === 'cloneInput' && isFastq) { loadedFastq=parseFastq(text); target.value=`FASTQ loaded: ${loadedFastq.length.toLocaleString()} reads\nClick Analyze construct to build a consensus.`; $('#cloneLength').textContent=`${loadedFastq.length.toLocaleString()} reads`; } else { loadedFastq=null; target.value=text; if(target.id==='cloneInput') $('#cloneLength').textContent=`${parseSequence(target.value).sequence.length.toLocaleString()} bp`; } }); }
function parseDelimitedRows(text, delimiter) {
  return String(text || '').trim().split(/\r?\n/).filter(Boolean).map(line => line.split(delimiter).map(cell => cell.trim().replace(/^"|"$/g, '')));
}
function insertRowsFromTable(rows) {
  const cleaned = rows.filter(row => row.some(cell => String(cell || '').trim()));
  if (!cleaned.length) return [];
  const headers = cleaned[0].map(cell => String(cell || '').trim().toLowerCase());
  const sequenceHeader = headers.findIndex(h => /^(sequence|seq|dna|insert|insert sequence)$/.test(h) || h.includes('sequence'));
  const nameHeader = headers.findIndex(h => /^(name|id|label|insert name)$/.test(h) || h.includes('name'));
  const hasHeader = sequenceHeader >= 0;
  const sequenceColumn = hasHeader ? sequenceHeader : cleaned[0].length > 1 ? 1 : 0;
  const nameColumn = hasHeader && nameHeader >= 0 ? nameHeader : cleaned[0].length > 1 ? 0 : -1;
  const data = cleaned.slice(hasHeader ? 1 : 0), insertsFromRows = [];
  for (let i = 0; i < data.length; i++) {
    const row = data[i], sequence = parseSequence(row[sequenceColumn] || '').sequence;
    if (!sequence) continue;
    insertsFromRows.push({ name: String(nameColumn >= 0 ? row[nameColumn] || `Insert ${i + 1}` : `Insert ${i + 1}`).trim(), sequence });
  }
  return insertsFromRows;
}
async function loadInsertTable(file) {
  const status = $('#insertFileStatus');
  try {
    let rows;
    if (/\.(csv|tsv)$/i.test(file.name)) {
      const text = await file.text();
      rows = parseDelimitedRows(text, /\.tsv$/i.test(file.name) ? '\t' : ',');
    } else {
      if (!window.XLSX) throw new Error('Excel support could not load. You can save the sheet as CSV and upload that instead.');
      const workbook = window.XLSX.read(await file.arrayBuffer(), { type: 'array' });
      const sheet = workbook.Sheets[workbook.SheetNames[0]];
      rows = window.XLSX.utils.sheet_to_json(sheet, { header: 1, raw: false, defval: '' });
    }
    const loaded = insertRowsFromTable(rows);
    if (!loaded.length) throw new Error('No DNA sequences were found. Use a Sequence column or put sequences in the second column.');
    inserts.innerHTML = '';
    loaded.forEach(item => addInsert(item.name, item.sequence));
    status.textContent = `${loaded.length} insert${loaded.length === 1 ? '' : 's'} loaded from ${file.name}. Uploading a new table replaces the current insert list.`;
  } catch (error) {
    status.textContent = error.message || 'The insert table could not be read.';
  }
}
loadFile($('#backboneFile'), $('#backboneInput'));
$('#insertFile').addEventListener('change', e => { const file=e.target.files[0]; if (file) loadInsertTable(file); });
document.querySelectorAll('.ghost').forEach(btn => btn.onclick = () => $(`#${btn.dataset.target}`).click());
$('#addInsert').onclick = () => addInsert(); $('#addClone').onclick = () => addClone(); $('#analyze').onclick = analyze;
$('#cloneBatchFile').addEventListener('change', async e => {
  const files = [...e.target.files]; if (!files.length) return;
  const current = [...document.querySelectorAll('.clone-sample')];
  if (current.length === 1 && !parseSequence(current[0].querySelector('.clone-input').value).sequence && !current[0]._reads) current[0].remove();
  for (const file of files) {
    const text = await file.text();
    const baseName = file.name.replace(/\.(fastq|fq|fasta|fa|txt)$/i, '') || `Clone ${document.querySelectorAll('.clone-sample').length + 1}`;
    const row = addClone(baseName);
    const textarea = row.querySelector('.clone-input'), length = row.querySelector('.clone-length');
    const isFastq = /\.(fastq|fq)$/i.test(file.name) || text.trimStart().startsWith('@') && text.split(/\r?\n/)[2] === '+';
    if (isFastq) { row._reads = parseFastq(text); textarea.value = `FASTQ loaded: ${row._reads.length} reads\nClick Check cloning result to build a consensus.`; length.textContent = `${row._reads.length.toLocaleString()} reads`; }
    else { textarea.value = text; length.textContent = `${parseSequence(text).sequence.length.toLocaleString()} bp`; }
  }
  e.target.value = '';
});
$('#demoButton').onclick = () => { loadedFastq=null; $('#backboneInput').value = '>backbone\nGGGAAACCCGGGTTTAAACCC'; $('#removedInput').value=''; inserts.innerHTML=''; addInsert('GFP insert','ATGAAACCCGGGTTT'); addInsert('Tag','GACTAG'); $('#clones').innerHTML=''; addClone('Clone 01','>clone-01\nGGGAAACCCGGGTTTAAACCCATGAAACCCGGGTTTGACTAG'); addClone('Clone 02','>clone-02\nGGGAAACCCGGGTTTAAACCCATGAAACCCGGGTTT'); };
$('#download').onclick = () => { if (!lastReport) return; const text = ['Plasmid Verify report', `Observed length: ${lastReport.observed.length} bp`, `Expected length: ${lastReport.construct.expected.length} bp`, `Differences: ${lastReport.diffs.length}`, '', ...lastReport.diffs.map(d => `${d.type}\tposition ${d.position}\t${d.expected} -> ${d.observed}`)].join('\n'); const a=document.createElement('a'); a.href=URL.createObjectURL(new Blob([text],{type:'text/plain'})); a.download='plasmid-verify-report.txt'; a.click(); URL.revokeObjectURL(a.href); };
addInsert('Insert 1');
addClone('Clone 01');
