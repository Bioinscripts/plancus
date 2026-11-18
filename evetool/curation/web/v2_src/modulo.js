
/* ================= curation_v2: the sequence viewer and the read heatmaps =================
   Sequences: one collapsible card per EVE locus, as the ViralQuest report shows its contigs. The
   header carries family, protein, length, verdict and the read evidence of every library mapped on
   the page; the body, read from the server when the card opens, shows the verdict, the genomic
   context with every read track, the domains, the reads library by library, the tblastn segments,
   the clusters and the sequence. Libraries: three heatmaps locus x library (DNA presence, RNA-seq
   transcription, small RNA with the piRNA signature) and the table of piRNA signatures. */
let ctxLarguraForcada = null;            // width for trilhasContexto when it draws inside a card
const svEstado = {abertos: new Set(), sel: new Set(), limite: 30, det: new Map(), lista: [], mapa: new Map(), iniciado: false};
const lbEstado = {rna: {}, srna: {}, pedindo: new Set(), iniciado: false, escolhida: false, carimbo: null};
const RNA_MIN_READS = 5, SRNA_MIN_READS = 10;
const PI_GLIFO = {pp: "●", prim: "▲", sec: "▼", si: "◆"};
const PI_ROT = {pp: "ping-pong", prim: "primary (5′U)", sec: "A10 only", si: "siRNA", size: "piRNA-sized, no bias",
                few: "few piRNA-sized reads", none: "no reads"};
const PI_ORDEM = ["pp", "prim", "sec", "si", "size", "few", "none"];
const lbCorta = (s, n) => { s = String(s ?? ""); return s.length > n ? s.slice(0, Math.max(1, n - 1)) + "…" : s; };

/* loci: the overlapping ORFs of one insertion together, through their representative ORF (as the server's
   loci_aprovados, so the heatmap rows and the matrices use the same key: species | locus) */
function svLoci(){
  const m = new Map();
  itens.forEach(i => { const k = (i.species || "") + "|" + (i.locus || i.query); if(!m.has(k)) m.set(k, []); m.get(k).push(i); });
  return [...m.entries()].map(([k, l]) => {
    const vivos = l.filter(i => !i.removed_flag), base = vivos.length ? vivos : l;
    const rep = base.find(i => i.locus_representative) || [...base].sort((a, b) => (+b.viral_bitscore || 0) - (+a.viral_bitscore || 0))[0];
    const span = vivos.length ? vivos : l;
    return {k, sp: rep.species || "", locus: rep.locus || rep.query, rep, orfs: l, queries: new Set(l.map(i => i.query)),
      start: Math.min(...span.map(i => i.start || 0)), end: Math.max(...span.map(i => i.end || 0)),
      rejeitado: !vivos.length, qualificado: vivos.some(i => i.qualified)};
  });
}

/* ---------- read evidence of a locus: DNA, RNA-seq, small RNA ---------- */
function lbIdx(m){ return m._idx || (m._idx = new Map((m.loci || []).map(x => [x.locus, x]))); }
function lbCel(m, L, nome){ const x = lbIdx(m).get(L.locus); return x ? (x.cells || {})[nome] || null : null; }
function lbDnaLibs(sp){
  const l = ((dnaDados && dnaDados.libraries) || {})[sp] || [];
  return [...l].sort((a, b) => (a.role === "assembly" ? 0 : 1) - (b.role === "assembly" ? 0 : 1)
    || (a.group || "").localeCompare(b.group || "") || a.name.localeCompare(b.name));
}
function lbDnaLinha(L){
  if(!dnaDados || !dnaDados.loci) return null;
  if(!dnaDados._idx){ dnaDados._idx = new Map(); dnaDados.loci.forEach(l => dnaDados._idx.set((l.species || "") + "|" + l.query, l)); }
  for(const q of L.queries){ const x = dnaDados._idx.get(L.sp + "|" + q); if(x) return x; }
  return dnaDados.loci.find(l => (l.species || "") === L.sp && l.contig === L.rep.contig && l.start < L.end && l.end > L.start) || null;
}
const LB_DNA = {pres: ["var(--s1)", 1], het: ["var(--s1)", .6], part: ["var(--s1)", .28], abs: ["var(--surface)", 1],
                exc: ["var(--s4)", 1], nd: ["var(--neutral-2)", 1]};
const LB_DNA_ROT = {pres: "present", het: "present · het", part: "partial", abs: "absent", exc: "excess", nd: "no data"};
function lbDnaClasse(c){ return c.c === "pres" && /het/.test(c.txt || "") ? "het" : c.c; }
function lbDnaCor(cl){ const [cor, op] = LB_DNA[cl] || LB_DNA.nd; return op >= 1 ? cor : `color-mix(in srgb, ${cor} ${Math.round(op * 100)}%, var(--surface))`; }
function lbDnaLegenda(){
  return `<div class="legend">${Object.keys(LB_DNA_ROT).map(k => `<span><i class="sw" style="background:${k === "nd"
    ? "repeating-linear-gradient(45deg,var(--neutral-2) 0 2px,var(--surface) 2px 4px)" : lbDnaCor(k)};border:1px solid var(--border)"></i>${LB_DNA_ROT[k]}</span>`).join("")}</div>`;
}
function lbDnaResumo(L){
  const libs = lbDnaLibs(L.sp);
  if(!libs.length) return null;
  const linha = lbDnaLinha(L), t = dnaLimiares();
  const ch = libs.map(b => ({b, c: dnaChamar(linha ? linha.dna[b.name] : null, b, t)}));
  return {linha, ch, n: libs.length, comDado: ch.filter(x => x.c.c !== "nd").length, pres: ch.filter(x => x.c.c === "pres").length};
}
function lbRnaResumo(L){
  const m = lbEstado.rna[L.sp];
  if(!m || !(m.samples || []).length) return null;
  const lst = m.samples.map(s => ({s, c: lbCel(m, L, s.name)}));
  return {lst, n: lst.length, transc: lst.filter(o => o.c && (o.c.reads || 0) >= RNA_MIN_READS).length,
    max: Math.max(0, ...lst.map(o => (o.c && o.c.mean_depth) || 0))};
}
function piLimiares(){
  const n = (id, d) => { const e = $(id), v = e ? parseFloat(e.value) : NaN; return isNaN(v) ? d : v; };
  return {min: n("#lb-pi-min", 20), u1: n("#lb-pi-u1", 0.5), a10: n("#lb-pi-a10", 0.4), z: n("#lb-pi-z", 2)};
}
/* piRNA signature of one locus in one library, from the piRNA-sized reads (24-31 nt) of each strand:
   primary = 5' U bias, secondary = A at position 10, ping-pong = z of the 10-nt 5'-5' overlap */
function piChamar(c, fitaEve, t){
  if(!c || !c.reads) return {c: "none", txt: "no small RNA read"};
  const p = c.pirna, si = c.cls === "siRNA";
  if(!p) return si ? {c: "si", txt: "siRNA (sRNAProfiler)"} : {c: "few", txt: "no piRNA-sized read"};
  const rel = f => fitaEve === "+" || fitaEve === "-" ? (f === fitaEve ? "sense" : "antisense") : (f === "+" ? "plus strand" : "minus strand");
  const s = ["+", "-"].map(f => ({f, rot: rel(f), n: (p[f] && p[f].reads) || 0, u: p[f] ? p[f].pct_1u : null, a: p[f] ? p[f].pct_10a : null}));
  const tot = s[0].n + s[1].n, base = {tot, s, z: p.pingpong_z, pares: p.pingpong_pairs10};
  if(tot < t.min) return {...base, c: si ? "si" : "few", txt: si ? "siRNA (sRNAProfiler)" : `${fmtN(tot, 0)} piRNA-sized reads (below ${t.min})`};
  const minF = Math.max(3, t.min / 4);
  const prim = s.filter(x => x.n >= minF && x.u != null && x.u >= t.u1);
  const sec = s.filter(x => x.n >= minF && x.a != null && x.a >= t.a10);
  const pp = (p.pingpong_z || 0) >= t.z && (p.pingpong_pairs10 || 0) > 0;
  const par = prim.find(a => sec.some(b => b.f !== a.f)), parB = par ? sec.find(b => b.f !== par.f) : null;
  if(pp) return {...base, c: "pp", txt: par ? `ping-pong cycle: primary (5′U) ${par.rot}, secondary (A10) ${parB.rot}, z ${p.pingpong_z}`
    : `ping-pong overlap (10 nt, z ${p.pingpong_z})${prim.length ? ` · 5′U ${prim.map(x => x.rot).join(" + ")}` : ""}${sec.length ? ` · A10 ${sec.map(x => x.rot).join(" + ")}` : ""}`};
  if(prim.length) return {...base, c: "prim", txt: `primary-like: 5′U ${prim.map(x => x.rot).join(" + ")}${sec.length ? ` · A10 ${sec.map(x => x.rot).join(" + ")}` : ""}, no ping-pong`};
  if(sec.length) return {...base, c: "sec", txt: `A10 ${sec.map(x => x.rot).join(" + ")}, no 5′U and no ping-pong`};
  if(si) return {...base, c: "si", txt: "siRNA (sRNAProfiler)"};
  return {...base, c: "size", txt: "piRNA-sized reads without 5′U or A10 bias"};
}
function lbSrnaResumo(L){
  const m = lbEstado.srna[L.sp];
  if(!m || !(m.libraries || []).length) return null;
  const t = piLimiares();
  const lst = m.libraries.map(b => { const c = lbCel(m, L, b.name); return {b, c, pi: piChamar(c, L.rep.strand, t)}; });
  const melhor = lst.map(o => o.pi.c).sort((a, b) => PI_ORDEM.indexOf(a) - PI_ORDEM.indexOf(b))[0];
  return {lst, n: lst.length, com: lst.filter(o => o.c && (o.c.reads || 0) >= SRNA_MIN_READS).length,
    reads: lst.reduce((s, o) => s + ((o.c && o.c.reads) || 0), 0), melhor};
}

/* ---------- loading: DNA measurements, RNA-seq and small RNA matrices, per species ---------- */
async function lbPedir(tipo, sp){
  const k = tipo + "|" + sp;
  if(lbEstado.pedindo.has(k)){ while(lbEstado.pedindo.has(k)) await new Promise(r => setTimeout(r, 150)); return; }
  lbEstado.pedindo.add(k);
  try{
    const r = await fetch(`api/${tipo === "rna" ? "rnaseq" : "srna"}-matrix?sp=` + encodeURIComponent(sp), {cache: "no-store"});
    const j = await r.json(); if(!r.ok || j.error) throw new Error(j.error || "HTTP " + r.status);
    lbEstado[tipo][sp] = j;
  }catch(e){ lbEstado[tipo][sp] = {error: String(e.message || e), loci: [], samples: [], libraries: []}; }
  finally{ lbEstado.pedindo.delete(k); }
}
async function lbDna(forcar){
  if(forcar && !dnaCarregando) dnaDados = null;
  if(!dnaDados && !dnaCarregando) await dnaCarregar();
  while(dnaCarregando) await new Promise(r => setTimeout(r, 150));
}
async function lbCarregarEspecie(sp, forcar){
  if(!vivo()) return;
  if(lbEstado.carimbo !== carimbo){ lbEstado.carimbo = carimbo; lbEstado.rna = {}; lbEstado.srna = {}; }   // curation changed: other loci
  const tarefas = [];
  if((rnaAmostras[sp] || []).length && (forcar || !lbEstado.rna[sp])) tarefas.push(lbPedir("rna", sp));
  if((srnaLibs[sp] || []).length && (forcar || !lbEstado.srna[sp])) tarefas.push(lbPedir("srna", sp));
  if(forcar || !dnaDados) tarefas.push(lbDna(forcar));
  await Promise.all(tarefas);
}

/* ---------- Sequences tab ---------- */
function svAbrir(){
  if(!svEstado.iniciado){ svEstado.iniciado = true; svLigar(); }
  svOpcoes(); svPintar();
  const sps = [...new Set(itens.map(i => i.species || ""))];
  Promise.all(sps.map(sp => lbCarregarEspecie(sp))).then(() => {
    if(document.body.dataset.tab !== "sequences") return;
    if($("#sv-reads").value || ["dna", "rna", "srna"].includes($("#sv-sort").value)) svPintar(); else svChipsAtualizar();
  });
}
function svOpcoes(){
  const enche = (sel, vals, todos, rot) => { const atual = sel.value;
    sel.innerHTML = `<option value="">${todos}</option>` + vals.map(v => `<option value="${escH(v)}">${escH(rot ? rot(v) : v)}</option>`).join("");
    if(vals.includes(atual)) sel.value = atual; };
  const sps = [...new Set(itens.map(i => i.species || ""))].sort();
  enche($("#sv-sp"), sps, "All species", sp => sp ? rotSpNome(sp) : "this genome");
  $("#sv-sp").hidden = sps.length < 2;
  enche($("#sv-fam"), [...new Set(itens.map(i => i.family).filter(f => f && f !== "n/d"))].sort(), "All families");
}
function svFiltrados(){
  const q = $("#sv-q").value.trim().toLowerCase(), sp = $("#sv-sp").value, fam = $("#sv-fam").value,
        est = $("#sv-state").value, rd = $("#sv-reads").value, num = v => +v || 0;
  let l = svLoci().filter(L => (!sp || L.sp === sp) && (!fam || L.rep.family === fam));
  l = l.filter(L => est === "all" ? true : est === "removed" ? L.rejeitado : est === "qualified" ? L.qualificado && !L.rejeitado
    : est === "undecided" ? !L.rejeitado && !L.qualificado : !L.rejeitado);
  if(q) l = l.filter(L => [L.rep.orf, L.sp, L.rep.family, L.rep.protein_class, L.rep.product, L.rep.contig, L.rep.ctx_cluster,
    L.rep.verdict, ...[...L.queries]].join(" ").toLowerCase().includes(q));
  if(rd) l = l.filter(L => rd === "dna" ? ((lbDnaResumo(L) || {}).pres || 0) > 0 : rd === "rna" ? ((lbRnaResumo(L) || {}).transc || 0) > 0
    : rd === "srna" ? ((lbSrnaResumo(L) || {}).com || 0) > 0 : (lbSrnaResumo(L) || {}).melhor === "pp");
  const ord = $("#sv-sort").value;
  if(ord === "pos") l.sort((a, b) => a.sp.localeCompare(b.sp) || String(a.rep.contig).localeCompare(String(b.rep.contig)) || a.start - b.start);
  else if(ord === "family") l.sort((a, b) => String(a.rep.family || "~").localeCompare(String(b.rep.family || "~")) || num(b.rep.viral_bitscore) - num(a.rep.viral_bitscore));
  else {
    const chave = {bits: L => num(L.rep.viral_bitscore), length: L => L.end - L.start, dna: L => (lbDnaResumo(L) || {}).pres || 0,
      rna: L => (lbRnaResumo(L) || {}).max || 0, srna: L => (lbSrnaResumo(L) || {}).reads || 0}[ord] || (L => num(L.rep.viral_bitscore));
    l.sort((a, b) => chave(b) - chave(a) || num(b.rep.viral_bitscore) - num(a.rep.viral_bitscore));
  }
  return l;
}
function svChips(L){
  const d = lbDnaResumo(L), r = lbRnaResumo(L), s = lbSrnaResumo(L), out = [];
  if(d) out.push(`<span class="chip${d.pres ? "" : " dim"}" title="DNA libraries where the insertion is present / libraries with data (thresholds of the Genomic DNA tab)"><span class="k">DNA</span>${d.pres}/${d.comDado}</span>`);
  if(r) out.push(`<span class="chip${r.transc ? "" : " dim"}" title="RNA-seq libraries with at least ${RNA_MIN_READS} reads inside the EVE, of ${r.n}; highest mean depth ${fmtN(r.max, 1)}×"><span class="k">RNA</span>${r.transc}/${r.n}</span>`);
  if(s) out.push(`<span class="chip${s.com ? "" : " dim"}" title="Small RNA reads over the EVE (1/NH weighted), ≥ ${SRNA_MIN_READS} in ${s.com} of ${s.n} libraries · strongest signature: ${PI_ROT[s.melhor]}"><span class="k">sRNA</span>${fmtN(s.reads, 0)}${PI_GLIFO[s.melhor] ? ` <span style="color:var(--s7)">${PI_GLIFO[s.melhor]}</span>` : ""}</span>`);
  return out.join(" ");
}
function svChipsAtualizar(){
  document.querySelectorAll("#sv-list .seq-card").forEach(c => { const L = svEstado.mapa.get(decodeURIComponent(c.dataset.k));
    const alvo = c.querySelector(".seq-reads"); if(L && alvo) alvo.innerHTML = svChips(L); });
}
function svCartao(L, maxBits){
  const i = L.rep, bits = +i.viral_bitscore || 0, len = L.end - L.start, gc = i.seq && i.seq.gc;
  const aberto = svEstado.abertos.has(L.k);
  return `<div class="seq-card${aberto ? " open" : ""}${L.rejeitado ? " rejected" : ""}" data-k="${encodeURIComponent(L.k)}">
    <div class="seq-head">
      <input type="checkbox" class="sv-sel"${svEstado.sel.has(L.k) ? " checked" : ""} title="Select for export" aria-label="Select ${escH(i.orf)}">
      <div class="seq-info">
      <span class="seq-id">${escH(i.orf)}</span>
      ${L.sp ? `<span class="seq-sp">${escH(rotSpNome(L.sp))}</span>` : ""}
      <span class="seq-meta"><b>Family</b>${escH(i.family && i.family !== "n/d" ? i.family : "—")}</span>
      <span class="seq-meta"><b>Protein</b>${escH(i.protein_class || "—")}</span>
      <span class="seq-meta">${len.toLocaleString()} nt</span>
      ${gc != null ? `<span class="seq-meta">GC ${(+gc).toFixed(1)}%</span>` : ""}
      ${L.orfs.length > 1 ? `<span class="seq-meta" title="ORFs of this locus (overlapping frames of one insertion)">${L.orfs.length} ORFs</span>` : ""}
      ${vdTag(i)}
      <span class="chip" title="Best viral bitscore (DIAMOND), against the highest of the list"><span class="k">BITS</span>${Math.round(bits)}<span class="meter"><i style="width:${(bits / maxBits * 100).toFixed(0)}%;background:var(--s1)"></i></span></span>
      <span class="seq-reads">${svChips(L)}</span>
      </div>
      <span class="seq-tail"><button class="mini sv-insp" title="Open the inspection (curation of this EVE) in a new tab">inspect</button><span class="seq-chev" aria-hidden="true">▶</span></span>
    </div>
    <div class="seq-body"${aberto ? "" : " hidden"}></div></div>`;
}
function svContagem(){
  const n = svEstado.lista.length, s = svEstado.sel.size;
  $("#sv-n").textContent = `${n.toLocaleString()} EVE loc${n === 1 ? "us" : "i"}${s ? ` · ${s} selected` : ""}`;
  $("#sv-export").disabled = !s;
}
function svPintar(){
  const el = $("#sv-list"); if(!el) return;
  const l = svFiltrados(), maxBits = Math.max(1, ...l.map(L => +L.rep.viral_bitscore || 0));
  svEstado.lista = l; svEstado.mapa = new Map(svLoci().map(L => [L.k, L]));
  const vis = l.slice(0, svEstado.limite);
  el.innerHTML = vis.length ? vis.map(L => svCartao(L, maxBits)).join("")
    : `<div class="card"><p class="sub" style="margin:0">${itens.length ? "No EVE matches the filters." : "Loading the candidates…"}</p></div>`;
  $("#sv-more").hidden = l.length <= vis.length;
  $("#sv-more").textContent = `show more (${(l.length - vis.length).toLocaleString()} left)`;
  svContagem();
  vis.forEach(L => { if(svEstado.abertos.has(L.k)) svCorpo(L); });
}
function svCartaoEl(k){ const e = encodeURIComponent(k); return [...document.querySelectorAll("#sv-list .seq-card")].find(c => c.dataset.k === e) || null; }
function svAlternar(k){
  const card = svCartaoEl(k), L = svEstado.mapa.get(k); if(!card || !L) return;
  const corpo = card.querySelector(".seq-body");
  if(svEstado.abertos.has(k)){ svEstado.abertos.delete(k); card.classList.remove("open"); corpo.hidden = true; return; }
  svEstado.abertos.add(k); card.classList.add("open"); corpo.hidden = false; svCorpo(L);
}
/* open a locus from elsewhere (a heatmap row, a cluster neighbour): its card, open, in view */
function svIr(k){
  const L = svLoci().find(x => x.k === k); if(!L) return;
  mostrarVista("sequences");
  if(!svFiltrados().some(x => x.k === k)){ $("#sv-q").value = L.rep.orf; $("#sv-sp").value = ""; $("#sv-fam").value = "";
    $("#sv-reads").value = ""; if(L.rejeitado) $("#sv-state").value = "all"; }
  const pos = svFiltrados().findIndex(x => x.k === k);
  if(pos >= svEstado.limite) svEstado.limite = pos + 1;
  svEstado.abertos.add(k); svPintar();
  const c = svCartaoEl(k); if(c) c.scrollIntoView({behavior: "smooth", block: "start"});
}
let svAtivos = 0;
async function svDetalhe(q){
  if(svEstado.det.has(q)) return svEstado.det.get(q);
  while(svAtivos >= 3) await new Promise(r => setTimeout(r, 120));       // a few at a time on "expand all"
  if(svEstado.det.has(q)) return svEstado.det.get(q);
  svAtivos++;
  try{
    const r = await fetch("api/detail?q=" + encodeURIComponent(q)); const d = await r.json();
    if(!r.ok || d.error) throw new Error(d.error || "HTTP " + r.status);
    svEstado.det.set(q, d); return d;
  }finally{ svAtivos--; }
}
const svSec = (titulo, ajuda) => `<div class="sec-label">${titulo}${ajuda ? `<span class="q" title="${escH(ajuda)}">?</span>` : ""}</div>`;
async function svCorpo(L){
  let card = svCartaoEl(L.k); if(!card) return;
  let corpo = card.querySelector(".seq-body");
  if(modo !== "live"){ corpo.innerHTML = `<p class="sub" style="margin:14px 0">Opening a card needs server.py (it reads the genome and the reads).</p>`; return; }
  if(!svEstado.det.has(L.rep.query)) corpo.innerHTML = `<p class="sub" style="margin:14px 0">Reading ${escH(L.rep.orf)}…</p>`;
  let d;
  try{ d = await svDetalhe(L.rep.query); }
  catch(e){ corpo.innerHTML = `<p class="sub" style="margin:14px 0">Could not read this EVE: ${escH(e.message || e)}</p>`; return; }
  card = svCartaoEl(L.k); if(!card || !svEstado.abertos.has(L.k)) return;       // redrawn or closed meanwhile
  corpo = card.querySelector(".seq-body");
  const W = Math.max(700, Math.floor(corpo.clientWidth - 40));
  corpo.innerHTML = svCorpoHTML(L, d, W);
  svLigarCorpo(corpo, L, d);
  pintarRnaseq(corpo.querySelector(".sv-ctx"));
}
function svCorpoHTML(L, d, W){
  const i = L.rep, nf = (x, k = 1) => x == null || x === "" || isNaN(+x) ? "—" : fmtN(+x, k);
  const met = [["viral bits", nf(i.viral_bitscore, 0)], ["cellular bits", nf(i.cellular_bitscore, 0)], ["identity", i.pident != null ? nf(i.pident, 1) + "%" : "—"],
    ["query cover", i.qcov != null ? nf(i.qcov, 0) + "%" : "—"], ["HMM score", nf(i.hmm_best_score || i.hmm_score, 1)],
    ["peptide", d.peptide_length ? d.peptide_length + " aa" : "—"], ["EVE GC", d.eve_gc != null ? d.eve_gc + "%" : "—"],
    ["strand", d.strand || i.strand || "—"]];
  let ctx;
  ctxLarguraForcada = W;
  try{ ctx = trilhasContexto(d); } finally { ctxLarguraForcada = null; }
  const fasta = `>${d.orf} ${d.contig}:${d.start + 1}-${d.end}(${d.strand})\n` + String(d.query_sequence || "").replace(/(.{70})/g, "$1\n");
  return `
    ${svSec("Verdict · evidence", "Explicit rules over the viral evidence and the evidence against it, as in the Inspection tab.")}
    ${vdCaixa(d.verdict || i)}
    <div class="met-row" style="margin-top:10px">${met.map(([k, x]) => `<div class="met"><div class="v">${escH(x)}</div><div class="k">${k}</div></div>`).join("")}</div>
    ${svSec("Genome map · context and reads", "The EVE in its genomic-context cluster: ORFs, other EVEs and host genes, and under them every RNA-seq and small RNA library mapped for this species (depth per strand). The same tracks as the inspection.")}
    <div class="sv-ctx">${ctx}</div>
    ${svSec("Domain architecture", "Every source of evidence placed on the ORF, in amino acids: viral alignment, viral HMM, CD-Search, Pfam and low complexity.")}
    <div>${arquiteturaDominios(d, W)}</div>
    ${svSec("Reads in the libraries", "DNA: presence of the insertion library by library. RNA-seq: reads inside the EVE. Small RNA: reads, class and piRNA signature. All libraries of the species, as in the Libraries tab.")}
    <div class="sec-grid">${svDnaBox(L)}${svRnaBox(L)}${svSrnaBox(L)}</div>
    ${svSec("tblastn segments", "Segments of the tblastn scan of the viral references (HSPs merged within the distance of the tblastn tab) that overlap this locus.")}
    ${svTbBox(L)}
    ${svSec("Clusters", "Genomic-context cluster: EVEs close together on the contig, possibly pieces of one viral genome. Cross-species cluster: EVEs similar between species (tBLASTx), possibly one insertion in a common ancestor.")}
    ${svClBox(L)}
    ${svSec("Sequence", "The ORF as read from the genome, on the strand of the ORF.")}
    <pre class="seq-fasta">${escH(fasta)}</pre>
    <div class="bar" style="margin-top:8px">
      <button class="mini sv-copy">Copy FASTA</button>
      ${botaoCopiar(d.species, d.contig, d.start, d.end, d.orf).replace(">copy<", ">copy genome region<")}
      <span style="flex:1"></span>
      <button class="mini primary sv-insp2">Open the inspection to curate</button></div>`;
}
function svLigarCorpo(corpo, L, d){
  const fasta = `>${d.orf}\n${d.query_sequence || ""}\n`;
  const c = corpo.querySelector(".sv-copy"); if(c) c.onclick = () => copiarTexto(fasta, c, "FASTA");
  const b = corpo.querySelector(".sv-insp2"); if(b) b.onclick = () => abrirInspecao(L.rep.query);
  corpo.querySelectorAll(".sv-goto").forEach(x => x.onclick = () => svIr(decodeURIComponent(x.dataset.k)));
}
function svDnaBox(L){
  const r = lbDnaResumo(L);
  if(!r) return `<div class="pbox"><h5>Genomic DNA</h5><p class="sub" style="margin:0">No DNA library for this species yet (Genomic DNA tab).</p></div>`;
  const grupos = new Map();
  r.ch.forEach(x => { const g = x.b.role === "assembly" ? "assembly reads" : (x.b.group || x.b.sample_type || "other");
    if(!grupos.has(g)) grupos.set(g, []); grupos.get(g).push(x); });
  return `<div class="pbox"><h5>Genomic DNA · present in ${r.pres} of ${r.comDado} with data</h5>${[...grupos.entries()].map(([g, l]) => {
    const pres = l.filter(x => x.c.c === "pres").length, com = l.filter(x => x.c.c !== "nd").length;
    return `<div class="lib-row"><span class="nm" title="${escH(g)}">${escH(g)}</span>
      <span>${l.map(x => { const cl = lbDnaClasse(x.c);
        return `<i class="libsq" style="background:${cl === "nd" ? "repeating-linear-gradient(45deg,var(--neutral-2) 0 2px,var(--surface) 2px 4px)" : lbDnaCor(cl)};border:1px solid var(--border)"
          title="${escH(x.b.name + ": " + x.c.txt)}"></i>`; }).join("")}</span>
      <span class="num">${pres}/${com}</span></div>`; }).join("")}${lbDnaLegenda()}</div>`;
}
function svRnaBox(L){
  const r = lbRnaResumo(L);
  if(!r) return `<div class="pbox"><h5>RNA-seq</h5><p class="sub" style="margin:0">${(rnaAmostras[L.sp] || []).length ? "Counting the reads…" : "No RNA-seq library for this species."}</p></div>`;
  const mx = Math.max(1e-9, r.max);
  return `<div class="pbox"><h5>RNA-seq · transcribed in ${r.transc} of ${r.n}</h5>${r.lst.map(o => { const c = o.c || {};
    return `<div class="lib-row"><span class="nm" title="${escH(o.s.name)} · library ${escH(o.s.library)}">${escH(o.s.name)}</span>
      <span class="barra" title="mean depth ${fmtN(c.mean_depth || 0, 1)}×"><i style="width:${((c.mean_depth || 0) / mx * 100).toFixed(1)}%;background:var(--s3)"></i></span>
      <span class="num" title="reads with aligned bases inside the EVE · mean depth · share of the EVE covered · depth on the + / − genome strand">${fmtN(c.reads || 0, 0)} reads · ${fmtN(c.mean_depth || 0, 1)}× · ${Math.round((c.breadth || 0) * 100)}%</span></div>`; }).join("")}</div>`;
}
function svSrnaBox(L){
  const r = lbSrnaResumo(L), m = lbEstado.srna[L.sp];
  if(!r) return `<div class="pbox"><h5>Small RNA</h5><p class="sub" style="margin:0">${(srnaLibs[L.sp] || []).length ? "Reading the small RNA libraries…" : "No small RNA library for this species."}</p></div>`;
  const melhor = r.lst.filter(o => o.c && o.c.reads).sort((a, b) => b.c.reads - a.c.reads)[0];
  const graf = melhor && melhor.c.sizes_plus ? `<p class="desc" style="margin:10px 0 4px">${escH(melhor.b.name)}, the library with most reads:
      read length by strand, and the 5′–5′ overlap of its piRNA-sized reads (10 nt = ping-pong)</p>
    <div class="srna-charts">${srnaTamanhos({...melhor.c, strand: L.rep.strand}, m)}${srnaSobreposicao({overlaps: (melhor.c.pirna || {}).overlaps || []})}</div>
    <p class="desc" style="margin:4px 0 0">${escH(melhor.pi.txt)}${melhor.pi.s ? " · " + melhor.pi.s.map(x => `${x.rot}: ${fmtN(x.n, 0)} piRNA-sized, 5′U ${x.u == null ? "—" : Math.round(x.u * 100) + "%"}, A10 ${x.a == null ? "—" : Math.round(x.a * 100) + "%"}`).join(" · ") : ""}</p>` : "";
  return `<div class="pbox"><h5>Small RNA · reads in ${r.com} of ${r.n}</h5>${r.lst.map(o => { const c = o.c || {};
    return `<div class="lib-row"><span class="nm" title="${escH(o.b.name)}">${escH(o.b.name)}</span>
      <span>${c.reads ? `${fmtN(c.reads, 0)} reads · ${fmtN(c.reads_rpm || 0, 2)} RPM · ${escH(c.cls || "—")}` : `<span class="desc">no read</span>`}</span>
      <span class="pi-tag ${o.pi.c}" title="${escH(o.pi.txt)}">${PI_GLIFO[o.pi.c] ? PI_GLIFO[o.pi.c] + " " : ""}${escH(PI_ROT[o.pi.c])}</span></div>`; }).join("")}${graf}</div>`;
}
function svTbBox(L){
  if(!varredura || !varredura.length) return `<p class="sub" style="margin:0">The tblastn scan has not been run (tblastn tab).</p>`;
  let segs = [];
  try{ segs = tbSegmentos(false).filter(s => (s.species || "") === L.sp && s.contig === L.rep.contig && s.start < L.end && s.end > L.start); }catch(e){}
  if(!segs.length) return `<p class="sub" style="margin:0">No tblastn segment overlaps this locus.</p>`;
  const cores = tbCores(segs);
  return segs.slice(0, 3).map(s => `<div class="pbox" style="margin-bottom:8px;background:var(--surface)"><h5>${escH(s.nome)} · ${escH(s.contig)}:${s.start.toLocaleString()}–${s.end.toLocaleString()} · ${s.hsps.length} HSP${s.hsps.length > 1 ? "s" : ""}</h5>${tbMapa(s, cores)}</div>`).join("")
    + (segs.length > 3 ? `<p class="desc">${segs.length - 3} more segment(s) in the tblastn tab.</p>` : "");
}
function svClBox(L){
  const partes = [], cx = L.rep.ctx_cluster;
  if(cx){
    const viz = svLoci().filter(o => o.sp === L.sp && o.k !== L.k && o.orfs.some(i => i.ctx_cluster === cx))
      .sort((a, b) => a.start - b.start);
    partes.push(`<div class="pbox"><h5>Genomic-context cluster · ${escH(cx)}</h5>${viz.length ? viz.map(o => `<div class="lib-row">
        <span class="nm"><code>${escH(o.rep.orf)}</code></span><span>${escH(o.rep.family || "—")} · ${escH(o.rep.protein_class || "")}</span>
        <span>${vdTag(o.rep)} <button class="mini sv-goto" data-k="${encodeURIComponent(o.k)}">open</button></span></div>`).join("")
      : `<p class="sub" style="margin:0">No other EVE in this cluster.</p>`}</div>`);
  } else partes.push(`<div class="pbox"><h5>Genomic-context cluster</h5><p class="sub" style="margin:0">Not in a context cluster (Clusters tab).</p></div>`);
  let xs = null;
  try{
    const um = (xsDados && xsDados.units_members) || {};
    xs = xsClusters().find(c => c.members.some(m => (um[m.eve] && um[m.eve].length ? um[m.eve] : [m.eve]).some(q => L.queries.has(q))));
  }catch(e){}
  if(xs) partes.push(`<div class="pbox"><h5>Cross-species cluster${xs.id != null ? " · " + escH(xs.id) : xs.name ? " · " + escH(xs.name) : ""}</h5>
      <p style="margin:0 0 6px;font-size:13px">${xs.sps.length} species · ${escH(xs.partilha)}${xs.idMax ? ` · best identity between species ${fmtN(xs.idMax, 0)}%` : ""}${
        xs.familias.length ? ` · ${escH(xs.familias.slice(0, 3).join(", "))}` : ""}</p>
      <div class="desc">${xs.sps.map(s => `<i>${escH(rotSpNome(s))}</i>`).join(", ")}</div></div>`);
  else partes.push(`<div class="pbox"><h5>Cross-species cluster</h5><p class="sub" style="margin:0">${xsDados ? "Not shared with another species." : "The cross-species clustering has not been run (Cross-species tab)."}</p></div>`);
  return `<div class="sec-grid">${partes.join("")}</div>`;
}
function svLigar(){
  $("#sv-list").addEventListener("click", ev => {
    const card = ev.target.closest(".seq-card"); if(!card) return;
    const k = decodeURIComponent(card.dataset.k), L = svEstado.mapa.get(k);
    const caixa = ev.target.closest(".sv-sel");
    if(caixa){ caixa.checked ? svEstado.sel.add(k) : svEstado.sel.delete(k); svContagem(); return; }
    if(ev.target.closest(".sv-insp")){ if(L) abrirInspecao(L.rep.query); return; }
    if(ev.target.closest(".seq-head")) svAlternar(k);
  });
  $("#sv-q").addEventListener("input", () => { svEstado.limite = 30; svPintar(); });
  ["#sv-sort", "#sv-sp", "#sv-fam", "#sv-state", "#sv-reads"].forEach(id => $(id).onchange = () => { svEstado.limite = 30; svPintar(); });
  $("#sv-more").onclick = () => { svEstado.limite += 30; svPintar(); };
  $("#sv-expand").onclick = () => { svEstado.lista.slice(0, svEstado.limite).forEach(L => svEstado.abertos.add(L.k)); svPintar(); };
  $("#sv-collapse").onclick = () => { svEstado.abertos.clear(); svPintar(); };
  $("#sv-selall").onclick = () => { svEstado.lista.forEach(L => svEstado.sel.add(L.k)); svPintar(); };
  $("#sv-desel").onclick = () => { svEstado.sel.clear(); svPintar(); };
  $("#sv-export").onclick = ev => { ev.stopPropagation(); $("#sv-export-menu").hidden = !$("#sv-export-menu").hidden; };
  document.addEventListener("click", ev => { if(!ev.target.closest(".menu")) $("#sv-export-menu").hidden = true; });
  $("#sv-export-menu").addEventListener("click", ev => { const b = ev.target.closest("button[data-x]"); if(!b) return;
    $("#sv-export-menu").hidden = true; svExportar(b.dataset.x); });
}
async function svExportar(x){
  const sel = svLoci().filter(L => svEstado.sel.has(L.k));
  if(!sel.length){ status("Select at least one EVE (the box at the left of each card).", "error"); return; }
  if(x === "gff"){          // the approved ORFs of each selected locus (all of them when the whole locus is rejected)
    orfsGFF(sel.flatMap(L => { const vivos = L.orfs.filter(i => !i.removed_flag); return vivos.length ? vivos : L.orfs; }), "eves_selected_orfs.gff3");
    return;
  }
  if(x === "tsv"){
    const linhas = [["species", "orf", "locus", "contig", "start", "end", "strand", "length", "family", "protein_class", "verdict",
      "viral_bitscore", "pident", "ctx_cluster", "dna_present", "dna_with_data", "rna_transcribed", "rna_libraries", "rna_max_mean_depth",
      "srna_reads", "srna_libraries_with_reads", "pirna_signature"].join("\t")];
    sel.forEach(L => { const d = lbDnaResumo(L), r = lbRnaResumo(L), s = lbSrnaResumo(L), i = L.rep;
      linhas.push([L.sp, i.orf, L.locus, i.contig, L.start + 1, L.end, i.strand, L.end - L.start, i.family, i.protein_class, i.verdict,
        i.viral_bitscore, i.pident, i.ctx_cluster, d ? d.pres : "", d ? d.comDado : "", r ? r.transc : "", r ? r.n : "",
        r ? r.max : "", s ? Math.round(s.reads * 100) / 100 : "", s ? s.com : "", s ? PI_ROT[s.melhor] : ""].map(v => v ?? "").join("\t")); });
    linhas.unshift(`# DNA calls with the thresholds of the Genomic DNA tab; RNA transcribed = ≥ ${RNA_MIN_READS} reads inside the EVE; piRNA calls: ${JSON.stringify(piLimiares())}`);
    baixarTSV(linhas, "eves_selected_reads.tsv"); status(`${sel.length} EVE loci exported.`, "ok"); return;
  }
  const regs = sel.map(L => ({sp: L.sp, c: L.rep.contig, s: L.start + 1, e: L.end, n: L.rep.orf,
    d: [L.rep.family, L.rep.protein_class].filter(v => v && v !== "n/d").join(" ")}));
  await fastaRegioes(regs, x === "fasta500" ? 500 : 0, "eves_selected.fasta", x !== "copy", null);
}

/* the "⋯" menu of the navigation bar: reload, theme, help, classic page (their own handlers stay) */
$("#b-more").onclick = ev => { ev.stopPropagation(); $("#nav-menu").hidden = !$("#nav-menu").hidden; };
$("#nav-menu").addEventListener("click", () => { $("#nav-menu").hidden = true; });
document.addEventListener("click", ev => { if(!ev.target.closest(".nav-menu")) $("#nav-menu").hidden = true; });
/* tabs that do not fit scroll sideways: fade the right edge while there is more to the right */
function navFade(){ const t = document.querySelector("header.top .tabs"); if(!t) return;
  t.classList.toggle("mais", t.scrollLeft + t.clientWidth < t.scrollWidth - 2); }
document.querySelector("header.top .tabs").addEventListener("scroll", navFade);
addEventListener("resize", navFade); requestAnimationFrame(navFade);

/* ---------- Libraries tab: heatmaps locus x library ---------- */
function lbAbrir(){
  if(!lbEstado.iniciado){ lbEstado.iniciado = true; lbLigar(); }
  const sps = [...new Set(itens.map(i => i.species || ""))].sort(), sel = $("#lb-sp"), atual = sel.value;
  const comLibs = sp => (rnaAmostras[sp] || []).length || (srnaLibs[sp] || []).length || (((dnaDados && dnaDados.libraries) || {})[sp] || []).length;
  sel.innerHTML = sps.map(sp => `<option value="${escH(sp)}">${escH(sp ? rotSpNome(sp) : "this genome")}${comLibs(sp) ? "" : " (no reads)"}</option>`).join("");
  if(sps.includes(atual) && (lbEstado.escolhida || comLibs(atual))) sel.value = atual;
  else { const com = sps.find(comLibs); if(com != null) sel.value = com; }
  sel.closest("label").hidden = sps.length < 2;
  lbPintar();
  const sp = sel.value;
  lbCarregarEspecie(sp).then(() => {
    if(!lbEstado.escolhida && !comLibs(sel.value)){ const com = sps.find(comLibs);
      if(com != null && com !== sel.value){ sel.value = com; return lbAbrir(); } }
    if(document.body.dataset.tab === "libraries") lbPintar();
  });
}
function lbLigar(){
  $("#lb-sp").onchange = () => { lbEstado.escolhida = true; lbPintar(); lbCarregarEspecie($("#lb-sp").value).then(lbPintar); };
  ["#lb-sort", "#lb-rna-v", "#lb-srna-v"].forEach(id => $(id).onchange = lbPintar);
  ["#lb-pi-min", "#lb-pi-u1", "#lb-pi-a10", "#lb-pi-z"].forEach(id => $(id).addEventListener("input", lbPintar));
  $("#lb-reload").onclick = async () => { const b = $("#lb-reload"); b.textContent = "…"; await rnaCarregarAmostras();
    await lbCarregarEspecie($("#lb-sp").value, true); b.textContent = "reload reads"; lbAbrir(); };
  $("#lb-dna-tsv").onclick = () => lbTSV("dna"); $("#lb-rna-tsv").onclick = () => lbTSV("rna"); $("#lb-srna-tsv").onclick = () => lbTSV("srna");
  $("#lb-pi").addEventListener("click", ev => { const tr = ev.target.closest("tr[data-k]"); if(tr) svIr(decodeURIComponent(tr.dataset.k)); });
}
function lbLinhas(){
  const sp = $("#lb-sp").value, ord = $("#lb-sort").value;
  const l = svLoci().filter(L => L.sp === sp && !L.rejeitado && L.rep.contig && L.end);
  if(ord === "family") l.sort((a, b) => String(a.rep.family || "~").localeCompare(String(b.rep.family || "~")) || a.start - b.start);
  else if(ord === "dna") l.sort((a, b) => ((lbDnaResumo(b) || {}).pres || 0) - ((lbDnaResumo(a) || {}).pres || 0) || a.start - b.start);
  else if(ord === "rna") l.sort((a, b) => ((lbRnaResumo(b) || {}).max || 0) - ((lbRnaResumo(a) || {}).max || 0) || a.start - b.start);
  else if(ord === "srna") l.sort((a, b) => ((lbSrnaResumo(b) || {}).reads || 0) - ((lbSrnaResumo(a) || {}).reads || 0) || a.start - b.start);
  else l.sort((a, b) => String(a.rep.contig).localeCompare(String(b.rep.contig)) || a.start - b.start);
  return l;
}
function lbPintar(){
  if(!$("#lb-sp").options.length) return;
  const sp = $("#lb-sp").value, loci = lbLinhas();
  $("#lb-n").textContent = `${loci.length.toLocaleString()} approved loci · ${sp ? rotSpNome(sp) : "this genome"}`;
  const nDna = lbDnaLibs(sp).length, mr = lbEstado.rna[sp], ms = lbEstado.srna[sp];
  const contar = f => loci.filter(f).length;
  const tiles = [[loci.length, "approved loci"], [nDna, "DNA libraries"]];
  if(nDna) tiles.push([contar(L => ((lbDnaResumo(L) || {}).pres || 0) > 0), "present in DNA"]);
  tiles.push([(rnaAmostras[sp] || []).length, "RNA-seq libraries"]);
  if(mr && mr.samples) tiles.push([contar(L => ((lbRnaResumo(L) || {}).transc || 0) > 0), "transcribed"]);
  tiles.push([(srnaLibs[sp] || []).length, "small RNA libraries"]);
  if(ms && ms.libraries) tiles.push([contar(L => ((lbSrnaResumo(L) || {}).com || 0) > 0), "with small RNA"],
                                    [contar(L => (lbSrnaResumo(L) || {}).melhor === "pp"), "piRNA ping-pong"]);
  $("#lb-tiles").innerHTML = tiles.map(([n, l]) => `<div class="tile"><div class="n">${(+n).toLocaleString()}</div><div class="l">${l}</div></div>`).join("");
  lbDnaHeat(loci); lbRnaHeat(loci); lbSrnaHeat(loci); lbPiTabela(loci);
}
/* one heatmap: rows = loci, columns = libraries (optionally grouped), a cell per pair */
function lbHeatmap(el, o){
  if(!o.colunas.length){ el.innerHTML = `<p class="sub" style="margin:0">${o.vazio}</p>`; return; }
  if(!o.linhas.length){ el.innerHTML = `<p class="sub" style="margin:0">No approved EVE in this genome.</p>`; return; }
  const larg = larguraCard(el, 600), nc = o.colunas.length;
  const rotW = Math.min(300, 16 + Math.max(...o.linhas.map(L => L.rep.orf.length)) * 6.3), famW = 128, extraW = o.extra ? 92 : 0;
  const cw = Math.max(12, Math.min(nc <= 12 ? 64 : 38, Math.floor((larg - rotW - famW - extraW - 20) / nc)));
  const ch = o.linhas.length > 80 ? 13 : o.linhas.length > 40 ? 15 : 18;
  const grupos = o.colunas.some(c => c.grupo);
  const topo = Math.min(150, 14 + Math.max(...o.colunas.map(c => Math.min(24, c.nome.length))) * 6.2 * 0.77);
  const y0 = (grupos ? 24 : 4) + topo, x0 = rotW + famW;
  const W = Math.max(larg, Math.ceil(x0 + nc * cw + extraW + 12)), H = Math.ceil(y0 + o.linhas.length * ch + 6);
  const id = "hatch-" + el.id;
  const p = [`<defs><pattern id="${id}" width="4" height="4" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
    <rect width="4" height="4" fill="var(--surface)"></rect><line x1="0" y1="0" x2="0" y2="4" stroke="var(--neutral-2)" stroke-width="1.6"></line></pattern></defs>`];
  if(grupos){
    let k = 0;
    while(k < nc){ const g = o.colunas[k].grupo || ""; let j = k; while(j + 1 < nc && (o.colunas[j + 1].grupo || "") === g) j++;
      const w = (j - k + 1) * cw;
      p.push(`<g><rect x="${x0 + k * cw + 1}" y="4" width="${w - 2}" height="15" rx="4" fill="var(--surface-2)" stroke="var(--border)"></rect>
        <text x="${x0 + k * cw + w / 2}" y="15" text-anchor="middle" font-size="10">${escH(lbCorta(g, Math.max(2, Math.floor(w / 6.2))))}</text>
        <title>${escH(g)}: ${j - k + 1} librar${j > k ? "ies" : "y"}</title></g>`);
      k = j + 1; }
  }
  o.colunas.forEach((c, k) => { const x = x0 + k * cw + cw / 2 + 3, y = y0 - 5;
    p.push(`<text x="${x}" y="${y}" transform="rotate(-50 ${x} ${y})" font-size="10.5"><title>${escH(c.nome)}${c.sub ? " · " + escH(c.sub) : ""}</title>${escH(lbCorta(c.nome, 24))}</text>`); });
  if(o.extra) p.push(`<text x="${x0 + nc * cw + 10}" y="${y0 - 6}" font-size="10" class="vlabel">${escH(o.extra.titulo)}</text>`);
  o.linhas.forEach((L, r) => {
    const y = y0 + r * ch, ym = y + ch / 2 + 3.5;
    const cel = o.colunas.map((c, k) => { const v = o.celula(L, k), x = x0 + k * cw;
      const fill = v.hachura ? `url(#${id})` : v.cor;
      return `<rect x="${x + 1}" y="${y + 1}" width="${cw - 2}" height="${ch - 2}" rx="2.5" fill="${fill}"${v.borda ? ' stroke="var(--border)"' : ""}><title>${escH(v.dica)}</title></rect>`
        + (v.glifo ? `<text x="${x + cw / 2}" y="${ym}" text-anchor="middle" font-size="${ch > 15 ? 10 : 8.5}" style="fill:${v.escuro ? "#fff" : "var(--text)"};pointer-events:none">${v.glifo}</text>` : ""); }).join("");
    p.push(`<g class="hm-row" data-k="${encodeURIComponent(L.k)}"><rect class="hm-bg" x="0" y="${y}" width="${x0 - 4}" height="${ch}" fill="transparent"></rect>
      <text class="rot" x="4" y="${ym}">${escH(L.rep.orf)}<title>${escH(L.rep.orf)} · ${escH(L.rep.contig)}:${(L.start + 1).toLocaleString()}–${L.end.toLocaleString()} · click to open its card</title></text>
      <text x="${rotW + 2}" y="${ym}" font-size="10.5" fill="var(--text-2)">${escH(lbCorta(L.rep.family && L.rep.family !== "n/d" ? L.rep.family : "—", 20))}</text>
      ${cel}${o.extra ? `<text x="${x0 + nc * cw + 10}" y="${ym}" font-size="10.5" class="vlabel">${escH(o.extra.valor(L))}</text>` : ""}</g>`);
  });
  el.innerHTML = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="${escH(o.titulo)}">${p.join("")}</svg>${o.legenda || ""}`;
  el.querySelectorAll("g.hm-row").forEach(g => g.addEventListener("click", () => svIr(decodeURIComponent(g.dataset.k))));
}
const lbVazioSp = sp => escH(sp ? rotSpNome(sp) : "this genome");
function lbGradiente(cor, mx, unidade, log){
  return `<div class="legend"><span>0</span><span class="grad" style="background:linear-gradient(90deg, color-mix(in srgb, ${cor} 12%, var(--surface)), ${cor})"></span>
    <span>${fmtN(mx, mx < 10 ? 2 : 0)} ${unidade}${log ? " (log scale)" : ""}</span>
    <span><i class="sw" style="background:var(--surface);border:1px solid var(--border)"></i>no read</span></div>`;
}
const lbTom = (cor, t) => `color-mix(in srgb, ${cor} ${Math.round(12 + 88 * Math.max(0, Math.min(1, t)))}%, var(--surface))`;
function lbDnaHeat(loci){
  const sp = $("#lb-sp").value, libs = lbDnaLibs(sp), t = dnaLimiares(), el = $("#lb-dna");
  if(!dnaDados){ el.innerHTML = `<p class="sub" style="margin:0">Reading the DNA libraries…</p>`; return; }
  const pop = libs.filter(b => b.role !== "assembly" && b.sample_type === "individual");
  lbHeatmap(el, {titulo: "Presence of each EVE in each DNA library", linhas: loci,
    colunas: libs.map(b => ({nome: b.name, grupo: b.role === "assembly" ? "assembly" : (b.group || b.sample_type || ""),
      sub: `${b.role === "assembly" ? "assembly reads" : b.sample_type}${b.genome_depth ? " · genome " + b.genome_depth + "×" : ""}`})),
    celula: (L, k) => { const linha = lbDnaLinha(L), m = linha ? linha.dna[libs[k].name] : null, c = dnaChamar(m, libs[k], t), cl = lbDnaClasse(c);
      return {cor: lbDnaCor(cl), borda: cl === "abs", hachura: cl === "nd",
        dica: `${L.rep.orf} · ${libs[k].name}: ${c.txt}\n` + (linha ? dnaDica(linha, m, libs[k], c) : "this locus is not among the DNA measurements")}; },
    extra: pop.length ? {titulo: "individuals", valor: L => { const r = lbDnaResumo(L); if(!r) return "—";
      const ind = r.ch.filter(x => pop.includes(x.b)), com = ind.filter(x => x.c.c !== "nd");
      return com.length ? `${ind.filter(x => x.c.c === "pres").length} / ${com.length}` : "—"; }} : null,
    legenda: lbDnaLegenda(),
    vazio: `No DNA library for ${lbVazioSp(sp)} yet: map reads in the <b>Genomic DNA</b> tab, or import the measurements made on a cluster.`});
}
function lbRnaHeat(loci){
  const sp = $("#lb-sp").value, m = lbEstado.rna[sp], campo = $("#lb-rna-v").value, el = $("#lb-rna");
  if(!(rnaAmostras[sp] || []).length){ el.innerHTML = `<p class="sub" style="margin:0">No RNA-seq library for ${lbVazioSp(sp)}: add a BAM or map reads from the read tracks of any EVE (a card or the inspection).</p>`; return; }
  if(!m){ el.innerHTML = `<p class="sub" style="margin:0">Counting the RNA-seq reads over every locus…</p>`; return; }
  if(m.error){ el.innerHTML = `<p class="sub" style="margin:0">Could not count the RNA-seq reads: ${escH(m.error)}</p>`; return; }
  const val = c => !c || c.error ? 0 : campo === "breadth" ? (c.breadth || 0) : (c[campo] || 0);
  const mx = Math.max(1e-9, ...loci.flatMap(L => m.samples.map(s => val(lbCel(m, L, s.name)))));
  const tr = v => campo === "breadth" ? v / mx : Math.log10(1 + v) / Math.log10(1 + mx);
  const un = {mean_depth: "× mean depth", reads: "reads", breadth: "of the EVE covered"}[campo];
  lbHeatmap(el, {titulo: "Transcription of each EVE in each RNA-seq library", linhas: loci,
    colunas: m.samples.map(s => ({nome: s.name, sub: "library " + s.library})),
    celula: (L, k) => { const s = m.samples[k], c = lbCel(m, L, s.name), v = val(c);
      if(c && c.error) return {cor: "var(--surface)", hachura: true, dica: `${s.name}: ${c.error}`};
      if(!v) return {cor: "var(--surface)", borda: true, dica: `${L.rep.orf} · ${s.name}: no read inside the EVE`};
      return {cor: lbTom("var(--s3)", tr(v)), dica: `${L.rep.orf} · ${s.name}\n${fmtN(c.reads || 0, 0)} reads inside the EVE · ${fmtN(c.spanning || 0, 0)} only spanning it by a splice gap
mean depth ${fmtN(c.mean_depth || 0, 2)}× (genome + strand ${fmtN(c.mean_plus || 0, 2)}, − strand ${fmtN(c.mean_minus || 0, 2)})
${Math.round((c.breadth || 0) * 100)}% of the EVE covered by at least one read`}; },
    extra: {titulo: "max", valor: L => { const v = Math.max(0, ...m.samples.map(s => val(lbCel(m, L, s.name))));
      return campo === "breadth" ? Math.round(v * 100) + "%" : fmtN(v, campo === "reads" ? 0 : 1) + (campo === "mean_depth" ? "×" : ""); }},
    legenda: lbGradiente("var(--s3)", campo === "breadth" ? 100 * mx : mx, campo === "breadth" ? "%" : un, campo !== "breadth")});
}
function lbSrnaDica(L, b, c, pi){
  const p = c.pirna || {}, s = pi.s || [];
  return [`${L.rep.orf} · ${b.name}`,
    `${fmtN(c.reads || 0, 1)} reads (1/NH weighted) · ${fmtN(c.reads_rpm || 0, 2)} RPM · ${c.distinct_seqs || 0} distinct sequences`,
    `sRNAProfiler class: ${c.cls || "—"}${c.confidence != null ? ` (confidence ${c.confidence})` : ""}${c.review ? " · review" : ""} · dominant length ${c.dicer_call || "—"} nt`,
    ...s.map(x => `${x.rot}: ${fmtN(x.n, 1)} piRNA-sized reads · 5′U ${x.u == null ? "—" : Math.round(x.u * 100) + "%"} · A10 ${x.a == null ? "—" : Math.round(x.a * 100) + "%"}`),
    p.pingpong_z != null ? `ping-pong (piRNA-sized): z ${p.pingpong_z}, ${fmtN(p.pingpong_pairs10 || 0, 1)} pairs overlapping 10 nt` : "",
    `call: ${pi.txt}`].filter(Boolean).join("\n");
}
const lbPiReads = c => c && c.pirna ? ((c.pirna["+"] || {}).reads || 0) + ((c.pirna["-"] || {}).reads || 0) : 0;
function lbSrnaHeat(loci){
  const sp = $("#lb-sp").value, m = lbEstado.srna[sp], campo = $("#lb-srna-v").value, el = $("#lb-srna"), t = piLimiares();
  if(!(srnaLibs[sp] || []).length){ el.innerHTML = `<p class="sub" style="margin:0">No small RNA library for ${lbVazioSp(sp)}: map trimmed reads or add a BAM from the read tracks of any EVE.</p>`; return; }
  if(!m){ el.innerHTML = `<p class="sub" style="margin:0">Reading the small RNA libraries over every locus…</p>`; return; }
  if(m.error){ el.innerHTML = `<p class="sub" style="margin:0">Could not read the small RNA libraries: ${escH(m.error)}</p>`; return; }
  const val = c => !c || c.error ? 0 : campo === "pirna" ? lbPiReads(c) : (c[campo] || 0);
  const mx = Math.max(1e-9, ...loci.flatMap(L => m.libraries.map(b => val(lbCel(m, L, b.name)))));
  const tr = v => Math.log10(1 + v) / Math.log10(1 + mx);
  lbHeatmap(el, {titulo: "Small RNA reads over each EVE in each library, with the piRNA signature", linhas: loci,
    colunas: m.libraries.map(b => ({nome: b.name, sub: b.lib_size ? `${Math.round(b.lib_size).toLocaleString()} mapped reads` : ""})),
    celula: (L, k) => { const b = m.libraries[k], c = lbCel(m, L, b.name), v = val(c);
      if(c && c.error) return {cor: "var(--surface)", hachura: true, dica: `${b.name}: ${c.error}`};
      if(!v) return {cor: "var(--surface)", borda: true, dica: `${L.rep.orf} · ${b.name}: no small RNA read`};
      const pi = piChamar(c, L.rep.strand, t), tt = tr(v);
      return {cor: lbTom("var(--s7)", tt), escuro: tt > .5, glifo: PI_GLIFO[pi.c] || "", dica: lbSrnaDica(L, b, c, pi)}; },
    extra: {titulo: "signature", valor: L => { const s = lbSrnaResumo(L); return s && PI_GLIFO[s.melhor] ? `${PI_GLIFO[s.melhor]} ${PI_ROT[s.melhor]}` : ""; }},
    legenda: lbGradiente("var(--s7)", mx, {reads_rpm: "RPM", reads: "reads", pirna: "piRNA-sized reads"}[campo], true)
      + `<div class="legend">${["pp", "prim", "sec", "si"].map(k => `<span><b style="color:var(--s7)">${PI_GLIFO[k]}</b>&nbsp;${PI_ROT[k]}</span>`).join("")}
         <span class="desc">thresholds above · hover a cell for the strands</span></div>`});
}
function lbSpark(v){
  const mx = Math.max(1e-9, ...v), W = 124, H = 26, bw = W / 20;
  return `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="5′–5′ overlap 1–20 nt">${v.map((c, k) => {
    const h = Math.max(c > 0 ? 1 : 0, c / mx * (H - 4));
    return `<rect x="${(k * bw + 0.5).toFixed(1)}" y="${(H - 2 - h).toFixed(1)}" width="${(bw - 1).toFixed(1)}" height="${h.toFixed(1)}"
      fill="${k === 9 ? "var(--s7)" : "var(--neutral)"}"><title>${k + 1} nt: ${fmtN(c, 1)} pairs</title></rect>`; }).join("")}
    <line x1="0" y1="${H - 1.5}" x2="${W}" y2="${H - 1.5}" stroke="var(--border)"></line></svg>`;
}
function lbPiTabela(loci){
  const sp = $("#lb-sp").value, m = lbEstado.srna[sp], t = piLimiares(), el = $("#lb-pi");
  if(!m || !m.libraries || m.error){ el.innerHTML = ""; return; }
  const linhas = [];
  loci.forEach(L => m.libraries.forEach(b => { const c = lbCel(m, L, b.name), pi = piChamar(c, L.rep.strand, t);
    if(pi.tot != null && pi.tot >= t.min) linhas.push({L, b, c, pi}); }));
  linhas.sort((a, b) => PI_ORDEM.indexOf(a.pi.c) - PI_ORDEM.indexOf(b.pi.c) || (b.pi.z || 0) - (a.pi.z || 0) || b.pi.tot - a.pi.tot);
  if(!linhas.length){ el.innerHTML = `<p class="sub" style="margin:10px">No locus with at least ${t.min} piRNA-sized reads in a library.</p>`; return; }
  const pc = v => v == null ? "—" : Math.round(v * 100) + "%";
  const lado = (s, fita) => { const a = s.find(x => x.f === (fita === "-" ? "-" : "+")), b = s.find(x => x !== a); return [a, b]; };
  const conhecida = f => f === "+" || f === "-";
  el.innerHTML = `<table class="pi-tab"><thead><tr><th>Locus</th><th>Viral family</th><th>Library</th><th>Signature</th>
      <th class="num" title="piRNA-sized reads (24–31 nt), 1/NH weighted, on the ORF's strand / the opposite strand">reads s / as</th>
      <th class="num" title="Share of the piRNA-sized reads of each strand starting with U">5′U s / as</th>
      <th class="num" title="Share of the piRNA-sized reads of each strand with an A at position 10">A10 s / as</th>
      <th class="num" title="z-score of the 10-nt 5′–5′ overlap between strands against 1–20 nt, piRNA-sized reads only">ping-pong z</th>
      <th>5′–5′ overlap, 1–20 nt</th></tr></thead><tbody>${linhas.map(({L, b, c, pi}) => {
      const [s, a] = lado(pi.s, L.rep.strand), nota = conhecida(L.rep.strand) ? "" : ' title="ORF strand unknown: + / − genome strand"';
      return `<tr data-k="${encodeURIComponent(L.k)}" style="cursor:pointer"><td><code>${escH(L.rep.orf)}</code></td><td>${escH(L.rep.family || "—")}</td>
        <td>${escH(b.name)}</td><td><span class="pi-tag ${pi.c}">${PI_GLIFO[pi.c] ? PI_GLIFO[pi.c] + " " : ""}${escH(PI_ROT[pi.c])}</span><br><span class="desc">${escH(pi.txt)}</span></td>
        <td class="num"${nota}>${fmtN(s.n, 0)} / ${fmtN(a.n, 0)}</td><td class="num">${pc(s.u)} / ${pc(a.u)}</td><td class="num">${pc(s.a)} / ${pc(a.a)}</td>
        <td class="num">${pi.z ?? "—"}<br><span class="desc">${fmtN(pi.pares || 0, 1)} pairs at 10</span></td>
        <td>${lbSpark(((c && c.pirna) || {}).overlaps || [])}</td></tr>`; }).join("")}</tbody></table>`;
}
function lbTSV(tipo){
  const sp = $("#lb-sp").value, loci = lbLinhas(), t = piLimiares(), saida = [];
  const cab = ["species", "orf", "locus", "contig", "start", "end", "strand", "family", "protein_class", "library"];
  if(tipo === "dna"){
    const libs = lbDnaLibs(sp), lt = dnaLimiares(); if(!libs.length) return status("No DNA library to export.", "error");
    saida.push([...cab, "group", "call", "eve_over_flank", "breadth", "eve_depth", "flank_depth", "span_left", "span_right", "empty_left", "empty_right"].join("\t"));
    loci.forEach(L => { const linha = lbDnaLinha(L); libs.forEach(b => { const m = linha ? linha.dna[b.name] : null, c = dnaChamar(m, b, lt);
      saida.push([sp, L.rep.orf, L.locus, L.rep.contig, L.start + 1, L.end, L.rep.strand, L.rep.family, L.rep.protein_class, b.name, b.group || b.role,
        c.txt, c.razao == null ? "" : c.razao.toFixed(3), m ? m.breadth : "", m ? m.eve_depth : "", m ? m.flank_depth : "",
        m ? m.span_left : "", m ? m.span_right : "", m ? m.empty_left : "", m ? m.empty_right : ""].map(v => v ?? "").join("\t")); }); });
  } else if(tipo === "rna"){
    const m = lbEstado.rna[sp]; if(!m || !(m.samples || []).length) return status("No RNA-seq counts to export.", "error");
    saida.push([...cab, "library_type", "reads_inside", "reads_spanning", "mean_depth", "breadth", "mean_plus", "mean_minus"].join("\t"));
    loci.forEach(L => m.samples.forEach(s => { const c = lbCel(m, L, s.name) || {};
      saida.push([sp, L.rep.orf, L.locus, L.rep.contig, L.start + 1, L.end, L.rep.strand, L.rep.family, L.rep.protein_class, s.name, s.library,
        c.reads ?? "", c.spanning ?? "", c.mean_depth ?? "", c.breadth ?? "", c.mean_plus ?? "", c.mean_minus ?? ""].join("\t")); }));
  } else {
    const m = lbEstado.srna[sp]; if(!m || !(m.libraries || []).length) return status("No small RNA reads to export.", "error");
    saida.push(`# piRNA calls: at least ${t.min} piRNA-sized reads; primary 5'U >= ${t.u1}; secondary A10 >= ${t.a10}; ping-pong z >= ${t.z}; strands relative to the ORF`);
    saida.push([...cab, "reads", "reads_rpm", "srnaprofiler_class", "confidence", "pirna_sense", "pirna_antisense", "pct_1u_sense", "pct_1u_antisense",
      "pct_10a_sense", "pct_10a_antisense", "pingpong_z_pirna", "pingpong_pairs10_pirna", "signature", "signature_detail"].join("\t"));
    loci.forEach(L => m.libraries.forEach(b => { const c = lbCel(m, L, b.name) || {}, pi = piChamar(c, L.rep.strand, t), p = c.pirna || {};
      const s = L.rep.strand === "-" ? "-" : "+", a = s === "+" ? "-" : "+", g = (f, k) => (p[f] || {})[k] ?? "";
      saida.push([sp, L.rep.orf, L.locus, L.rep.contig, L.start + 1, L.end, L.rep.strand, L.rep.family, L.rep.protein_class, b.name,
        c.reads ?? 0, c.reads_rpm ?? 0, c.cls ?? "", c.confidence ?? "", g(s, "reads"), g(a, "reads"), g(s, "pct_1u"), g(a, "pct_1u"),
        g(s, "pct_10a"), g(a, "pct_10a"), p.pingpong_z ?? "", p.pingpong_pairs10 ?? "", PI_ROT[pi.c], pi.txt].join("\t")); }));
  }
  baixarTSV(saida, `eve_${tipo}_${(sp || "genome").replace(/\W+/g, "_")}.tsv`);
}
