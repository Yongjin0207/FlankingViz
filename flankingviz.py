#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FlankingViz: local flanking-sequence mapping and junction-aware visualization.

Existing profile and setup.sh formats are preserved. Optional full-vector FASTA
or GenBank improves boundary evidence; vector-free reads remain supported with
labelled BLAST-based boundary estimates or interval-based annotations.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ----------------------------
# helpers
# ----------------------------

DNA_RE = re.compile(r"[^ACGTNacgtn]+")

def read_any_sequence(path: Path) -> str:
    """
    Accepts:
      - raw .seq (no header)
      - FASTA (with '>')
    Returns uppercase A/C/G/T/N only.
    """
    return read_sequence(path)

def write_fasta(tmp_fa: Path, name: str, seq: str) -> None:
    with tmp_fa.open("w", encoding="utf-8") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i+70] + "\n")

def run_cmd(cmd: List[str]) -> Tuple[int, str, str]:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stdout, p.stderr


# ----------------------------
# GFF gene model parsing
# ----------------------------

@dataclass
class Gene:
    chr: str
    start: int
    end: int
    strand: str
    display_id: str   # Os01t... or LOC_Os... etc (what we show)
    gene_id: str      # Os01g... or LOC...
    source: str       # RAP_ANN, RAP_PRED, MSU
    note: str = ""

    @property
    def size_bp(self) -> int:
        return int(self.end - self.start + 1)

def parse_attr(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in attr.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out

def load_genes_rap(path: Path, source: str) -> Dict[str, List[Gene]]:
    """
    RAP locus.gff style:
      chr01 irgsp1_locus gene 2983 10815 . + . ID=Os01g0100100; ... Transcript variants=Os01t0100100-01
    We show transcript variant (Os..t..) in display_id to match RAPDB UI.
    """
    genes: Dict[str, List[Gene]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attr = cols[:9]
            if ftype != "gene":
                continue
            try:
                st = int(start); en = int(end)
            except ValueError:
                continue
            a = parse_attr(attr)
            gid = a.get("ID", a.get("Name", "")).strip()
            tv = a.get("Transcript variants", "").split(",")[0].strip()
            # Some have "Os01t....-00" etc
            disp = tv if tv else gid
            note = a.get("Note", "").replace("%2C", ",")
            g = Gene(chr=seqid, start=st, end=en, strand=strand if strand in "+-" else "+",
                     display_id=disp, gene_id=gid, source=source, note=note)
            genes.setdefault(seqid, []).append(g)
    # sort for each chr
    for k in genes:
        genes[k].sort(key=lambda x: x.start)
    return genes

def load_genes_msu(path: Path) -> Dict[str, List[Gene]]:
    """
    MSU/UGA gene-only gff3 (already converted to chr01.. and only gene lines):
      chr01 MSU_osa1r7 gene 2903 10817 . + . ID=LOC_Os01g01010;Name=LOC_Os01g01010;...
    display_id = LOC_Os.. (same)
    """
    genes: Dict[str, List[Gene]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attr = cols[:9]
            if ftype != "gene":
                continue
            try:
                st = int(start); en = int(end)
            except ValueError:
                continue
            a = parse_attr(attr)
            gid = a.get("ID", a.get("Name", "")).strip()
            note = a.get("Note", "").replace("%20", " ")
            g = Gene(chr=seqid, start=st, end=en, strand=strand if strand in "+-" else "+",
                     display_id=gid, gene_id=gid, source="MSU", note=note)
            genes.setdefault(seqid, []).append(g)
    for k in genes:
        genes[k].sort(key=lambda x: x.start)
    return genes




def load_genes_generic_gff3(
    path: Path,
    source: str,
    feature: str = "gene",
    attr_priority: Optional[List[str]] = None,
) -> Dict[str, List[Gene]]:
    """
    Generic GFF3 gene parser for non-rice references (e.g., Soybean, Maize).

    - Uses rows where type == `feature` (default: "gene").
    - display_id is chosen by attribute priority (default: ["Name","gene","ID"]).
    - gene_id is set to ID if present, else display_id.
    """
    if attr_priority is None:
        attr_priority = ["Name", "gene", "ID"]

    genes: Dict[str, List[Gene]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attr = cols[:9]
            if ftype != feature:
                continue
            try:
                st = int(start); en = int(end)
            except ValueError:
                continue
            a = parse_attr(attr)

            disp = ""
            for k in attr_priority:
                v = a.get(k, "").strip()
                if v:
                    disp = v
                    break
            if not disp:
                disp = a.get("ID", "").strip() or a.get("Name", "").strip() or ""

            gid = a.get("ID", "").strip() or disp

            note = ""
            if "Note" in a:
                note = a.get("Note", "").replace("%2C", ",").replace("%20", " ")
            g = Gene(
                chr=seqid,
                start=st,
                end=en,
                strand=strand if strand in "+-" else "+",
                display_id=disp if disp else gid,
                gene_id=gid,
                source=source,
                note=note,
            )
            genes.setdefault(seqid, []).append(g)

    for k in genes:
        genes[k].sort(key=lambda x: x.start)
    return genes


def load_profile_yaml(path: Path) -> Dict[str, object]:
    """
    Lightweight YAML loader:
    - Prefers PyYAML if available.
    - Falls back to a minimal parser for the simple key/value + list format used in FlankingViz profiles.
    """
    try:
        import yaml  # type: ignore
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        if not isinstance(data, dict):
            raise ValueError("Profile YAML must be a mapping (dict).")
        return data
    except Exception:
        # Minimal YAML for:
        # key: value
        # key:
        #   - item1
        #   - item2
        out: Dict[str, object] = {}
        cur_key: Optional[str] = None
        cur_list: Optional[List[str]] = None

        for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("- "):
                if cur_key is None:
                    continue
                if cur_list is None:
                    cur_list = []
                    out[cur_key] = cur_list
                cur_list.append(line[2:].strip().strip("'\""))
                continue

            # new key
            if ":" in line:
                k, v = line.split(":", 1)
                k = k.strip()
                v = v.strip()
                cur_key = k
                cur_list = None
                if v == "":
                    out[k] = []
                    cur_list = out[k]  # type: ignore
                else:
                    out[k] = v.strip("'\"")
        return out
# ----------------------------
# Nearest gene selection
# ----------------------------

def overlaps(g: Gene, pos: int) -> bool:
    # legacy point-overlap
    return g.start <= pos <= g.end

def overlaps_interval(g: Gene, a: int, b: int) -> bool:
    """Return True if gene overlaps [a,b] (inclusive)."""
    return not (g.end < a or g.start > b)

def nearest_left(genes: List[Gene], pos: int) -> Optional[Tuple[Gene, int]]:
    best = None
    for g in genes:
        if g.end <= pos:
            dist = pos - g.end
            if best is None or dist < best[1]:
                best = (g, dist)
    return best

def nearest_right(genes: List[Gene], pos: int) -> Optional[Tuple[Gene, int]]:
    best = None
    for g in genes:
        if g.start >= pos:
            dist = g.start - pos
            if best is None or dist < best[1]:
                best = (g, dist)
    return best

def pick_overlap(genes: List[Gene], pos: int) -> Tuple[List[Gene], Optional[Gene]]:
    ovs = [g for g in genes if overlaps(g, pos)]
    if not ovs:
        return [], None
    # pick representative overlap gene (closest boundary)
    def score(g: Gene) -> int:
        return min(abs(pos - g.start), abs(g.end - pos))
    rep = sorted(ovs, key=lambda g: (score(g), g.size_bp))[0]
    return ovs, rep

def nearest_left_interval(genes: List[Gene], a: int) -> Optional[Tuple[Gene, int]]:
    """Nearest gene fully on the left of interval: gene.end < a; dist = a - gene.end"""
    best = None
    for g in genes:
        if g.end < a:
            dist = a - g.end
            if best is None or dist < best[1]:
                best = (g, dist)
    return best

def nearest_right_interval(genes: List[Gene], b: int) -> Optional[Tuple[Gene, int]]:
    """Nearest gene fully on the right of interval: gene.start > b; dist = gene.start - b"""
    best = None
    for g in genes:
        if g.start > b:
            dist = g.start - b
            if best is None or dist < best[1]:
                best = (g, dist)
    return best

def pick_overlap_interval(genes: List[Gene], a: int, b: int) -> Tuple[List[Gene], Optional[Gene]]:
    ovs = [g for g in genes if overlaps_interval(g, a, b)]
    if not ovs:
        return [], None
    # representative: closest boundary to interval edges
    def score(g: Gene) -> int:
        return min(abs(a - g.start), abs(a - g.end), abs(b - g.start), abs(b - g.end))
    rep = sorted(ovs, key=lambda g: (score(g), g.size_bp))[0]
    return ovs, rep


def union_genes(gene_maps: List[Dict[str, List[Gene]]], chr_id: str) -> List[Gene]:
    out: List[Gene] = []
    for m in gene_maps:
        out.extend(m.get(chr_id, []))
    # No need to sort fully for nearest scan; but helps determinism
    out.sort(key=lambda g: (g.start, g.end, g.display_id))
    return out


# ----------------------------
# BLASTN (RAPDB-like options)
# ----------------------------

@dataclass
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"

def parse_blast6(text: str) -> List[BlastHit]:
    hits: List[BlastHit] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        try:
            hit = BlastHit(
                qseqid=cols[0],
                sseqid=cols[1],
                pident=float(cols[2]),
                length=int(cols[3]),
                qstart=int(cols[6]),
                qend=int(cols[7]),
                sstart=int(cols[8]),
                send=int(cols[9]),
                evalue=float(cols[10]),
                bitscore=float(cols[11]),
            )
            hits.append(hit)
        except Exception:
            continue
    # Best first: higher bitscore, lower evalue
    hits.sort(key=lambda h: (-h.bitscore, h.evalue, -h.length))
    return hits

def is_noncanonical_rice_seqid(seqid: str) -> bool:
    """
    Exclude non-main rice chromosome hits only.
    Keep chr01..chr12, filter out records like ChrSy/Syng_*.
    """
    s = (seqid or "").strip().lower()
    if not s:
        return False
    if re.fullmatch(r"chr(0[1-9]|1[0-2])", s):
        return False
    if s == "chrsy" or s.startswith("syng_"):
        return True
    return False

def blast_best_hit(db: str, seq_name: str, seq: str, tmp_dir: Path) -> Optional[BlastHit]:
    qfa = tmp_dir / f"{seq_name}.fa"
    write_fasta(qfa, seq_name, seq)

    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    cmd = [
        "blastn",
        "-task", "blastn",
        "-db", db,
        "-query", str(qfa),
        "-max_target_seqs", "100",
        "-evalue", "0.1",
        "-dust", "yes",
        "-soft_masking", "true",
        "-outfmt", outfmt,
    ]
    rc, stdout, stderr = run_cmd(cmd)
    if rc != 0:
        raise RuntimeError(f"blastn failed for {seq_name}\nCMD: {' '.join(cmd)}\nSTDERR:\n{stderr}")
    hits = parse_blast6(stdout)

    # Rice-specific conservative filter: ignore non-canonical ChrSy / Syng_* hits
    db_name = Path(db).name.lower()
    if hits and db_name == "irgsp1.0":
        canonical_hits = [h for h in hits if not is_noncanonical_rice_seqid(h.sseqid)]
        if canonical_hits:
            return canonical_hits[0]

    return hits[0] if hits else None

def infer_insertion_pos(hit: BlastHit) -> int:
    """Legacy helper for an explicit query-start assumption.
    On both strands, sstart is the subject base corresponding to qstart.
    Automatic junction evidence is evaluated by infer_junction in this file.
    """
    # Both strands: sstart corresponds to qstart. This is an explicit side assumption.
    return hit.sstart


# ----------------------------
# HTML generation (single file)
# ----------------------------



def inject_ui_patch_js(template_html: str) -> str:
    """UI patch (no dependency on template internals)
    - Double-click list entries to rename (stored in localStorage, reapplied periodically)
    - Export CSV button downloads embedded CSV
    """
    patch = r'''<script>
(function(){
  const LS_KEY = "flankingviz_name_map_v1";

  function loadMap(){ try{ return JSON.parse(localStorage.getItem(LS_KEY)||"{}")||{}; }catch(e){ return {}; } }
  function saveMap(m){ try{ localStorage.setItem(LS_KEY, JSON.stringify(m||{})); }catch(e){} }

  function findListItems(){
    // Prefer explicit list containers if present
    const roots = [
      document.querySelector("#list"),
      document.querySelector("#recordList"),
      document.querySelector(".sidebar"),
      document.querySelector(".left"),
      document
    ].filter(Boolean);

    let best = [];
    for(const root of roots){
      const cand = Array.from(root.querySelectorAll("li, .item, .list-item, .row, .record, .entry, div"))
        .filter(el=>{
          const t = (el.textContent||"").trim();
          if(!t) return false;
          if(t.length > 140) return false; // avoid huge containers
          const tag = (el.tagName||"").toLowerCase();
          if(tag==="script"||tag==="style") return false;
          // heuristic: looks like a row
          const cls = (el.className||"").toString().toLowerCase();
          if(cls.includes("item")||cls.includes("row")||cls.includes("list")||cls.includes("record")||cls.includes("entry")) return true;
          if(el.getAttribute("role")==="button") return true;
          if(root && root.id==="list") return true;
          return false;
        });
      if(cand.length > best.length) best = cand;
      if(root && root.id==="list" && cand.length) break;
    }
    return best;
  }

  function applyMap(){
    const m = loadMap();
    const items = findListItems();
    for(const el of items){
      const raw = (el.getAttribute("data-raw-title") || el.textContent || "").trim();
      if(!raw) continue;
      if(!el.getAttribute("data-raw-title")) el.setAttribute("data-raw-title", raw);
      const mapped = m[raw];
      if(mapped && mapped.trim().length){
        el.textContent = mapped;
      }
    }
  }

  function renameEl(el){
    if(!el) return;
    const raw = (el.getAttribute("data-raw-title") || el.textContent || "").trim();
    if(!raw) return;
    if(!el.getAttribute("data-raw-title")) el.setAttribute("data-raw-title", raw);

    const m = loadMap();
    const cur = m[raw] || raw;
    const nn = prompt("표시 이름을 수정하세요", cur);
    if(nn === null) return;
    const name = (nn||"").trim();
    if(!name){
      delete m[raw];
      el.textContent = raw;
    }else{
      m[raw] = name;
      el.textContent = name;
    }
    saveMap(m);
  }

  function addExportBtn(){
    if(window.__FV_EXPORT_CSV__) return;
    window.__FV_EXPORT_CSV__ = true;

    const btn = document.createElement("button");
    btn.textContent = "Export CSV";
    btn.title = "Download the embedded CSV";
    btn.style.position = "fixed";
    btn.style.right = "14px";
    btn.style.bottom = "14px";
    btn.style.zIndex = "9999";
    btn.style.padding = "10px 12px";
    btn.style.borderRadius = "12px";
    btn.style.border = "1px solid rgba(0,0,0,0.15)";
    btn.style.background = "white";
    btn.style.boxShadow = "0 6px 16px rgba(0,0,0,0.12)";
    btn.style.fontSize = "13px";
    btn.style.cursor = "pointer";

    btn.addEventListener("click", function(){
      try{
        const csv = window.__EMBEDDED_CSV__;
        if(!csv || !csv.length){ alert("No embedded CSV found."); return; }
        const blob = new Blob([csv], {type:"text/csv;charset=utf-8"});
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = "flankingviz.csv";
        document.body.appendChild(a);
        a.click();
        a.remove();
        setTimeout(()=>URL.revokeObjectURL(url), 800);
      }catch(e){
        alert("Export failed.");
      }
    });
    document.body.appendChild(btn);
  }

  function bind(){
    applyMap();
    setInterval(applyMap, 800);

    // Capture dblclick so existing handlers won't swallow it
    function looksLikeMainTitle(el){
      try{
        if(!el) return false;
        const t = (el.textContent||"").trim();
        if(!t) return false;
        if(t.length > 160) return false;
        // Heuristic pattern used by this viewer's main title
        if(t.includes(" | line ") && (t.includes("Chr ") || t.includes("Chr\t") || t.includes("Chr")) ) return true;
        if(t.startsWith("seqs | line ") && t.includes("| Chr")) return true;
        return false;
      }catch(e){ return false; }
    }

    function renameTitleEl(el){
      if(!el) return;
      const raw = (el.getAttribute("data-raw-title") || el.textContent || "").trim();
      if(!raw) return;
      if(!el.getAttribute("data-raw-title")) el.setAttribute("data-raw-title", raw);
      const m = loadMap();
      const cur = m[raw] || raw;
      const nn = prompt("Please edit the title", cur);
      if(nn===null) return;
      const name = (nn||"").trim();
      if(!name){
        delete m[raw];
        saveMap(m);
        el.textContent = raw;
        return;
      }
      m[raw] = name;
      saveMap(m);
      el.textContent = name;
    }

    // Capture dblclick so existing handlers won't swallow it
    document.addEventListener("dblclick", function(ev){
      try{
        // 1) Main title (works for both HTML text and SVG <text>)
        let node = ev.target;
        while(node && node !== document.body){
          if(looksLikeMainTitle(node)){
            renameTitleEl(node);
            return;
          }
          node = node.parentElement;
        }

        // 2) List items
        const items = findListItems();
        node = ev.target;
        while(node && node !== document.body){
          if(items.includes(node)){
            renameEl(node);
            return;
          }
          node = node.parentElement;
        }
      }catch(e){}
    }, true);

    addExportBtn();
  }

  if(document.readyState === "loading"){
    document.addEventListener("DOMContentLoaded", bind);
  }else{
    bind();
  }
})();
</script>'''
    if "</body>" in template_html:
        return template_html.replace("</body>", patch + "\n</body>")
    return template_html + patch


def load_template(template_path: Path) -> str:
    return template_path.read_text(encoding="utf-8", errors="ignore")

def embed_csv_into_html(template_html: str, csv_text: str) -> str:
    # Safe JS string literal: escape backslashes, backticks, and ${ to avoid template literal injection
    safe = csv_text.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")
    inject = f"<script>\nwindow.__EMBEDDED_CSV__ = `{safe}`;\n</script>\n"
    # insert right after <body> for early availability
    return template_html.replace("<body>", "<body>\n" + inject, 1)


# ----------------------------
# Main
# ----------------------------

# ----------------------------
# Junction and candidate analysis
# ----------------------------

import csv
import hashlib
import io
import json
import os
import platform
import re
import subprocess
import tempfile
import time
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from pathlib import Path

ENGINE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
LEGACY = ['Species','line','Chr','insertion_strand','left_id','left_strand',
          'left_dist_bp','left_size_bp','right_id','right_strand','right_dist_bp','right_size_bp','Note']
EXTRA = ['candidate_id','candidate_rank','candidate_count','reported_candidate_count','mapping_status','junction_status',
         'reference_seqid',
         'insertion_pos','junction_min','junction_max','coordinate_basis','classification',
         'alignment_strand','alignment_start','alignment_end','query_start','query_end',
         'identity_pct','alignment_length','bitscore','evalue','vector_id','vector_query_start',
         'vector_query_end','junction_query_side','query_gap_bp','overlap_gene_ids',
         'overlap_gene_sources','alignment_overlap_gene_ids','left_gene_id','left_source',
         'left_gene_start','left_gene_end','left_relation','left_status','left_dist_min',
         'left_dist_max','right_gene_id','right_source','right_gene_start','right_gene_end',
         'right_relation','right_status','right_dist_min','right_dist_max','feature_context',
         'warnings','pair_status','annotation_basis','annotation_start','annotation_end',
         'boundary_inference_reason','query_length','unaligned_query_start_bp','unaligned_query_end_bp']
FIELDS = LEGACY + EXTRA
DNA = frozenset('ACGTRYSWKMBDHVN')

def read_sequence(path):
    text = Path(path).read_text(encoding='utf-8-sig').strip()
    if not text:
        return ''
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    headers = [i for i, line in enumerate(lines) if line.startswith('>')]
    if headers and (len(headers) != 1 or headers[0] != 0):
        raise ValueError(f'{path}: one sequence per input file is required; multi-record FASTA is not concatenated.')
    sequence = re.sub(r'\s+', '', ''.join(lines[1:] if headers else lines)).upper()
    invalid = sorted(set(sequence) - DNA)
    if invalid:
        raise ValueError(f'{path}: invalid sequence characters {invalid}; no bases were silently removed.')
    return sequence

def read_vectors(path):
    """FASTA (one/many records) or GenBank ORIGIN; CM5 must be exported first."""
    path = Path(path)
    if path.suffix.lower()=='.cm5':
        raise ValueError('Export the Clone Manager CM5 vector as FASTA or GenBank before registration.')
    text = path.read_text(encoding='utf-8-sig')
    records = []
    if text.lstrip().startswith('LOCUS'):
        for block in text.split('//'):
            if not block.strip():
                continue
            locus = re.search(r'^LOCUS\s+(\S+)', block, re.M)
            origin = re.search(r'^ORIGIN\s*\n(.*)', block, re.M | re.S)
            if not locus or not origin:
                raise ValueError('GenBank vector requires LOCUS and ORIGIN sequence.')
            records.append((locus[1], re.sub(r'[\s\d]', '', origin[1]).upper()))
    else:
        for block in text.split('>'):
            if not block.strip():
                continue
            if not text.lstrip().startswith('>'):
                raise ValueError('Vector file must be FASTA or GenBank. Export CM5 as FASTA/GenBank.')
            lines = block.splitlines()
            records.append((lines[0].strip().split()[0], re.sub(r'\s+', '', ''.join(lines[1:])).upper()))
    if not records or len({n for n, _ in records}) != len(records):
        raise ValueError('Vector records require nonempty, unique identifiers.')
    for name, seq in records:
        if len(seq) < 20 or set(seq) - DNA:
            raise ValueError(f'Invalid vector sequence: {name}')
    return records

@dataclass
class Hit:
    query: str
    chrom: str
    identity: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    score: float
    qseq: str = ''
    sseq: str = ''
    @property
    def strand(self):
        return '+' if self.sstart <= self.send else '-'
    @property
    def start(self):
        return min(self.sstart, self.send)
    @property
    def end(self):
        return max(self.sstart, self.send)

def blast_batch(records, db, tmp, args, vector=False):
    result = {name: [] for name, _ in records}
    for offset in range(0, len(records), args.batch_size):
        part = records[offset:offset+args.batch_size]
        query = tmp / ('vector_query.fa' if vector else 'host_query.fa')
        query.write_text(''.join(f'>{name}\n{seq}\n' for name,seq in part), encoding='ascii')
        command = ['blastn','-task','blastn-short' if vector else 'blastn','-db',str(db),
                   '-query',str(query),'-max_target_seqs',str(args.max_hits),'-evalue','0.1',
                   '-max_hsps',str(args.max_hsps),
                   '-dust','no' if vector else 'yes','-soft_masking','true',
                   '-num_threads',str(args.threads),'-outfmt',
                   '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq']
        proc = subprocess.run(command, capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(f'BLAST failed: {proc.stderr.strip()}')
        for line in proc.stdout.splitlines():
            c = line.split('\t')
            if len(c) != 14 or c[0] not in result:
                raise RuntimeError('Unexpected BLAST output; analysis stopped.')
            result[c[0]].append(Hit(c[0],c[1],float(c[2]),int(c[3]),int(c[6]),int(c[7]),
                                     int(c[8]),int(c[9]),float(c[10]),float(c[11]),c[12],c[13]))
    for hits in result.values():
        hits.sort(key=lambda h:(-h.score,h.evalue,-h.length,h.chrom,h.start,h.qstart,h.qend))
    return result

def distinct_loci(hits):
    """Collapse redundant HSPs only when both genomic and query spans overlap 80%."""
    selected = []; bins={}
    for h in hits:
        redundant = False
        nearby=set()
        keys=[(h.chrom,h.strand,b) for b in range(h.start//1000,h.end//1000+1)]
        for key in keys: nearby.update(bins.get(key,[]))
        for index in sorted(nearby):
            prior=selected[index]
            host_overlap = max(0,min(h.end,prior.end)-max(h.start,prior.start)+1)
            query_overlap = max(0,min(h.qend,prior.qend)-max(h.qstart,prior.qstart)+1)
            if (host_overlap >= .8*min(h.end-h.start+1,prior.end-prior.start+1)
                    and query_overlap >= .8*min(h.qend-h.qstart+1,prior.qend-prior.qstart+1)):
                redundant = True
                break
        if not redundant:
            for key in keys: bins.setdefault(key,[]).append(len(selected))
            selected.append(h)
    return selected

def mapped_host_bases(hit, query_low, query_high):
    if not hit.qseq:
        return [hit.sstart+(q-hit.qstart)*(1 if hit.strand=='+' else -1)
                for q in range(query_low,query_high+1)]
    query=hit.qstart;subject=hit.sstart;step=1 if hit.strand=='+' else -1;positions=[]
    for q,s in zip(hit.qseq,hit.sseq):
        if q!='-' and s!='-' and query_low<=query<=query_high:positions.append(subject)
        if q!='-':query+=1
        if s!='-':subject+=step
    return positions

def infer_unaligned_tail(hit, sequence, host_hits, args):
    """Heuristic query-end evidence only; does not identify the tail as vector DNA."""
    if sequence is None:
        return None, 'sequence_extent_unavailable'
    minimum=getattr(args,'min_unaligned_tail',20)
    terminal=getattr(args,'max_terminal_overhang',10)
    before,after=hit.qstart-1,len(sequence)-hit.qend
    if before>=minimum and after<=terminal:
        side,lo,hi,pos='start',1,hit.qstart-1,hit.sstart
    elif after>=minimum and before<=terminal:
        side,lo,hi,pos='end',hit.qend+1,len(sequence),hit.send
    else:
        return None, 'both_ends_unaligned' if before>terminal and after>terminal else 'no_clear_unaligned_tail'
    tail=sequence[lo-1:hi]
    if sum(c in 'ACGT' for c in tail)<minimum or sum(c not in 'ACGT' for c in tail)>len(tail)*.1:
        return None, 'ambiguous_bases_in_tail'
    # Another host alignment in this tail contradicts a simple non-host -> host transition.
    for other in host_hits:
        if other.identity>=args.min_identity and other.length>=args.min_host_anchor:
            overlap=max(0,min(hi,other.qend)-max(lo,other.qstart)+1)
            if overlap>=minimum:
                return None, 'additional_host_alignment_in_tail'
    return (side,pos), 'one_sided_unaligned_tail'


def infer_junction(hit, vector_hits, args, sequence=None, host_hits=()):
    result = dict(junction_status='unresolved_no_vector' if not args.vector else 'unresolved_vector_not_found',
                  insertion_pos='',junction_min=hit.start,junction_max=hit.end,
                  coordinate_basis='alignment_span_context_only',vector_id='',vector_query_start='',
                  vector_query_end='',junction_query_side='',query_gap_bp='',
                  boundary_inference_reason='',query_length=len(sequence) if sequence is not None else '',
                  unaligned_query_start_bp=hit.qstart-1,
                  unaligned_query_end_bp=len(sequence)-hit.qend if sequence is not None else '')
    # Require independent, adjacent alignment blocks. Overlap and filler remain uncertain.
    choices = []
    for v in vector_hits:
        if v.length < args.min_vector_anchor or v.identity < args.min_vector_identity:
            continue
        if v.qstart < hit.qstart and v.qend < hit.qend:
            gap, side, pos = hit.qstart-v.qend-1, 'start', hit.sstart
        elif v.qend > hit.qend and v.qstart > hit.qstart:
            gap, side, pos = v.qstart-hit.qend-1, 'end', hit.send
        else:
            continue
        if -args.max_junction_overlap <= gap <= args.max_junction_gap:
            choices.append((abs(gap),-v.score,side,pos,gap,v))
    if choices:
        choices.sort(key=lambda x:(x[0],x[1],x[2],x[5].chrom))
        _,_,side,pos,gap,v = choices[0]
        competing = {c[2] for c in choices if c[0] == choices[0][0] and -c[1] >= v.score*args.ambiguity_ratio}
        result.update(vector_id=v.chrom,vector_query_start=v.qstart,vector_query_end=v.qend,
                      junction_query_side=side,query_gap_bp=gap)
        if len(competing) > 1:
            result['junction_status'] = 'unresolved_multiple_junction_sides'
        elif gap != 0:
            result['junction_status'] = 'unresolved_query_gap' if gap>0 else 'unresolved_query_overlap'
            if gap<0:
                low,high=(hit.qstart,min(hit.qend,v.qend+1)) if side=='start' else (max(hit.qstart,v.qstart-1),hit.qend)
                positions=mapped_host_bases(hit,low,high)
                if positions:
                    result.update(junction_min=min(positions),junction_max=max(positions),
                                  coordinate_basis='vector_host_overlap_candidate_interval')
        else:
            result.update(junction_status='vector_supported',insertion_pos=pos,junction_min=pos,junction_max=pos,
                          coordinate_basis='1_based_adjacent_aligned_host_base')
    elif args.junction_side != 'auto':
        pos = hit.sstart if args.junction_side == 'start' else hit.send
        result.update(junction_status='assumed_query_'+args.junction_side,insertion_pos=pos,junction_min=pos,
                      junction_max=pos,junction_query_side=args.junction_side,
                      coordinate_basis='1_based_host_base_explicit_side_assumption')
    elif sequence is not None:
        inferred,reason=infer_unaligned_tail(hit,sequence,host_hits,args)
        result['boundary_inference_reason']=reason
        if inferred:
            side,pos=inferred
            result.update(junction_status='estimated_query_'+side,insertion_pos=pos,
                          junction_min=pos,junction_max=pos,junction_query_side=side,
                          coordinate_basis='1_based_host_base_blast_tail_estimate')
    return result


def gene_priority(g):
    # Display-only tie break. All overlaps are retained; priority does not establish biology.
    return ({'RAP_ANN':0,'MSU':1,'RAP_PRED':2}.get(g.source,3),g.gene_id,g.start,g.end)

class GeneIndex:
    """Coordinate index; preserves every annotation model and tie-breaking rule."""
    def __init__(self, genes):
        self.by_start=sorted(genes,key=lambda g:(g.start,g.end,gene_priority(g)))
        self.by_end=sorted(genes,key=lambda g:(g.end,gene_priority(g)))
        self.starts=[g.start for g in self.by_start];self.ends=[g.end for g in self.by_end]
        self.prefix_end=[];largest=0
        for g in self.by_start:
            largest=max(largest,g.end);self.prefix_end.append(largest)
    def overlapping(self,a,b):
        i=bisect_right(self.starts,b)-1;result=[]
        while i>=0 and self.prefix_end[i]>=a:
            if self.by_start[i].end>=a:result.append(self.by_start[i])
            i-=1
        return sorted(result,key=gene_priority)
    def left(self,a,legacy=False):
        i=bisect_left(self.ends,a)-1
        if i<0:return None
        j=bisect_left(self.ends,self.ends[i])
        key=(lambda g:(g.start,g.end,g.display_id)) if legacy else gene_priority
        return min(self.by_end[j:i+1],key=key)
    def right(self,b,legacy=False):
        i=bisect_right(self.starts,b)
        if i>=len(self.starts):return None
        j=bisect_right(self.starts,self.starts[i])
        key=(lambda g:(g.start,g.end,g.display_id)) if legacy else gene_priority
        return min(self.by_start[i:j],key=key)

def annotate(genes, hit, junction):
    index=genes if isinstance(genes,GeneIndex) else GeneIndex(genes)
    point = junction['insertion_pos'] != ''
    a,b = (junction['insertion_pos'],junction['insertion_pos']) if point else (hit.start,hit.end)
    overlaps = index.overlapping(a,b)
    status = 'intragenic' if overlaps else 'intergenic'
    left = index.left(a,legacy=not point)
    right = index.right(b,legacy=not point)
    # Interval fallback reproduces the original representative-gene rule.
    if overlaps:
        left = overlaps[0] if point else min(overlaps,key=lambda g:(
            min(abs(a-g.start),abs(a-g.end),abs(b-g.start),abs(b-g.end)),
            g.size_bp,g.start,g.end,g.display_id))
    basis = ('blast_estimated_junction' if junction['junction_status'].startswith('estimated_') else
             'vector_supported_junction' if junction['junction_status']=='vector_supported' else
             'explicit_side_assumption') if point else 'alignment_span_fallback'
    result = dict(classification=status,annotation_basis=basis,annotation_start=a,annotation_end=b,overlap_gene_ids=';'.join(g.gene_id for g in overlaps),
                  overlap_gene_sources=';'.join(g.source for g in overlaps),
                  alignment_overlap_gene_ids=';'.join(g.gene_id for g in index.overlapping(hit.start,hit.end)),
                  feature_context='gene_body_subfeature_unavailable' if overlaps else 'not_assigned')
    for side,g in [('left',left),('right',right)]:
        for field in ['id','strand','dist_bp','size_bp','gene_id','source','gene_start','gene_end','relation','dist_min','dist_max']:
            result[f'{side}_{field}'] = ''
        result[f'{side}_status'] = 'no_gene_on_side' if g is None else 'annotated_gene'
        if g is None:
            continue
        within = g.start<=b and g.end>=a
        low,high = (0,0) if within else ((a-g.end,b-g.end) if side=='left' else (g.start-b,g.start-a))
        relation = ('within_gene' if point else 'overlaps_gene') if within else (
            ('downstream' if g.strand=='+' else 'upstream') if side=='left' else
            ('upstream' if g.strand=='+' else 'downstream'))
        result.update({f'{side}_id':g.display_id,f'{side}_gene_id':g.gene_id,f'{side}_source':g.source,
                       f'{side}_strand':g.strand,f'{side}_size_bp':g.size_bp,f'{side}_gene_start':g.start,
                       f'{side}_gene_end':g.end,f'{side}_dist_bp':('X' if within else low),
                       f'{side}_dist_min':low,f'{side}_dist_max':high,f'{side}_relation':relation})
    return result

def make_rows(name, species, hits, vectors, genes_by_chr, args, rice=False, sequence=None):
    warnings=[]
    all_host_hits=hits
    by_subject={}
    for h in hits: by_subject[h.chrom]=by_subject.get(h.chrom,0)+1
    if any(n>=args.max_hsps for n in by_subject.values()):warnings.append('hsp_limit_reached')
    if rice:
        primary=[h for h in hits if not (h.chrom.lower()=='chrsy' or h.chrom.lower().startswith('syng_'))]
        if primary:
            if len(primary)!=len(hits): warnings.append('nonprimary_rice_hits_excluded')
            hits=primary
    if len({h.chrom for h in hits}) >= args.max_hits:
        warnings.append('target_limit_reached')
    qualified=[h for h in hits if h.length>=args.min_host_anchor and h.identity>=args.min_identity]
    eligible=qualified or hits[:1]
    candidates=distinct_loci(eligible)
    if not candidates:
        row=dict.fromkeys(FIELDS,'')
        row.update(Species=species,line=name,Chr='-',candidate_id='1',candidate_rank=1,
            candidate_count=0,mapping_status='no_hit',junction_status='unresolved',classification='unresolved',Note='No BLAST hit')
        return [row]
    if not qualified: warnings.append('below_host_alignment_threshold')
    close=sum(h.score>=candidates[0].score*args.ambiguity_ratio for h in candidates)
    mapping='low_quality' if not qualified else ('ambiguous' if close>1 else 'unique_best')
    count=len(candidates)
    if count>args.max_candidates:warnings.append('candidate_output_limit_reached')
    if mapping=='unique_best' and any('limit_reached' in w for w in warnings):mapping='search_limited'
    candidates=candidates[:args.max_candidates]
    rows=[]
    for i,h in enumerate(candidates):
        j=infer_junction(h,vectors,args,sequence,all_host_hits)
        if not qualified:
            j.update(boundary_inference_reason='below_host_alignment_threshold',junction_status='unresolved_low_quality',insertion_pos='',junction_min=h.start,junction_max=h.end,
                     coordinate_basis='alignment_span_context_only')
        row=dict.fromkeys(FIELDS,'')
        chromosome=(h.chrom[3:].lstrip('0') or h.chrom) if h.chrom.startswith('chr') else h.chrom
        row.update(Species=species,line=name,Chr=chromosome,reference_seqid=h.chrom,candidate_id=str(i+1),candidate_rank=i+1,
                   candidate_count=count,reported_candidate_count=len(candidates),mapping_status=mapping if i==0 else (
                       'alternative_near_best' if h.score>=candidates[0].score*args.ambiguity_ratio else 'alternative_lower_score'),
                   insertion_strand=h.strand,alignment_strand=h.strand,alignment_start=h.start,alignment_end=h.end,
                   query_start=h.qstart,query_end=h.qend,identity_pct=h.identity,alignment_length=h.length,
                   bitscore=h.score,evalue=h.evalue,warnings=';'.join(warnings),**j)
        row.update(annotate(genes_by_chr.get(h.chrom,[]),h,j))
        row['Note']=f"Candidate {i+1}/{len(candidates)}; {row['mapping_status']}; {j['junction_status']}; {row['annotation_basis']}; alignment {h.chrom}:{h.start}-{h.end}({h.strand})"
        rows.append(row)
    return rows

def csv_text(rows, fields=FIELDS):
    output=io.StringIO(newline='')
    writer=csv.DictWriter(output,fieldnames=fields,extrasaction='ignore',lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return output.getvalue()

def check_pairs(path, rows, max_distance):
    by_name={r['line']:r for r in rows}; results=[];seen=set()
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        if not {'pair_id','line_a','line_b'}.issubset(reader.fieldnames or []):
            raise ValueError('Pair CSV requires pair_id,line_a,line_b columns.')
        for p in reader:
            if not p['pair_id'] or p['pair_id'] in seen: raise ValueError('Pair IDs must be unique and nonempty.')
            seen.add(p['pair_id'])
            if p['line_a']==p['line_b'] or any(p[k] not in by_name for k in ['line_a','line_b']):
                raise ValueError(f"Pair {p['pair_id']}: use two distinct existing sequence names.")
            a,b=by_name[p['line_a']],by_name[p['line_b']];distance=''
            if any(r['mapping_status']!='unique_best' for r in [a,b]): status='unresolved_mapping'
            elif a.get('reference_seqid',a['Chr'])!=b.get('reference_seqid',b['Chr']): status='discordant_chromosomes'
            elif any(r['junction_status']!='vector_supported' for r in [a,b]): status='unresolved_junction'
            else:
                distance=abs(a['insertion_pos']-b['insertion_pos'])
                status='nearby_same_chromosome' if distance<=max_distance else 'distant_same_chromosome'
            results.append(dict(pair_id=p['pair_id'],line_a=p['line_a'],line_b=p['line_b'],
                status=status,distance_bp=distance,screening_distance_bp=max_distance,
                interpretation='Coordinate screening only; does not establish insert architecture or copy number.'))
            for r in [a,b]: r['pair_status']=';'.join(filter(None,[r['pair_status'],p['pair_id']+':'+status]))
    return results

def safe_json(value):
    return json.dumps(value,ensure_ascii=False).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')

def reference_provenance(profile, api):
    if not profile: return {'status':'legacy_paths_in_parameters'}
    profile=Path(profile);cfg=api.load_profile_yaml(profile);root=profile.parent.parent
    result={'profile_sha256':hashlib.sha256(profile.read_bytes()).hexdigest(),'annotations':[]}
    files=cfg.get('gff_files',[])
    if isinstance(files,str):files=[files]
    for rel in files:
        path=Path(rel)
        if not path.is_absolute():path=root/path
        h=hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
        result['annotations'].append(dict(path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest()))
    manifest=root/'validation/reference_manifest.json'
    if manifest.is_file():
        data=json.loads(manifest.read_text())
        result['reference_inventory']=data.get('profiles',{}).get(str(cfg.get('name','')), {})
        result['inventory_checked_utc']=data.get('checked_utc','')
        result['inventory_note']='Inventory is from the stated validation time. Current annotation hashes are recorded separately.'
    return result

def run_pipeline(args, gene_maps, species, db, profile, api):
    start=time.perf_counter()
    out=Path(args.out).resolve();out.parent.mkdir(parents=True,exist_ok=True)
    if out.suffix.lower() != '.html':
        raise ValueError('--out must end in .html')
    paths={suffix:out.with_suffix(suffix) for suffix in ['.html','.csv','.candidates.csv','.json','.run.json']}
    if args.pairs: paths['.pairs.csv']=out.with_suffix('.pairs.csv')
    if not args.overwrite and any(p.exists() for p in paths.values()):
        raise SystemExit('Output already exists. Use a new --out name or explicitly pass --overwrite.')
    lock=out.with_suffix('.lock')
    try: lock_fd=os.open(str(lock),os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    except FileExistsError: raise SystemExit(f'Output is in use (lock {lock}). Choose a different output path.')
    os.close(lock_fd)
    try:
        files=sorted(p for p in Path(args.seqs_dir).resolve().iterdir() if p.is_file() and p.suffix.lower() in {'.seq','.fa','.fasta','.txt'})
        if not files: raise ValueError('No sequence files found.')
        if len({p.stem for p in files})!=len(files): raise ValueError('Duplicate input names: file stems must be unique.')
        if any(any(ord(c)<32 for c in p.stem) for p in files):raise ValueError('Input names must not contain control characters.')
        records=[];info=[];short=[]
        for i,p in enumerate(files):
            sequence=read_sequence(p);qid=f'q{i:07d}'
            info.append(dict(query_id=qid,line=p.stem,file=p.name,length=len(sequence),sha256=hashlib.sha256(sequence.encode()).hexdigest()))
            if len(sequence)<args.min_len: short.append(p.stem)
            else: records.append((qid,sequence))
        by_chr={}
        for gm in gene_maps:
            for chrom,genes in gm.items(): by_chr.setdefault(chrom,[]).extend(genes)
        for chrom in by_chr:
            # Exact duplicates within one source only; distinct annotation models remain visible.
            by_chr[chrom]=GeneIndex(list({(g.source,g.gene_id,g.start,g.end,g.strand):g for g in by_chr[chrom]}.values()))
        with tempfile.TemporaryDirectory(prefix='.flankingviz-',dir=out.parent) as folder:
            tmp=Path(folder)
            hits=blast_batch(records,db,tmp,args) if records else {}
            vector_hits={};vector_meta=[]
            if args.vector:
                vectors=read_vectors(args.vector);vf=tmp/'vector.fa'
                vf.write_text(''.join(f'>v{i:05d}\n{s}\n' for i,(_,s) in enumerate(vectors)),encoding='ascii')
                vector_names={f'v{i:05d}':n for i,(n,_) in enumerate(vectors)}
                proc=subprocess.run(['makeblastdb','-in',str(vf),'-dbtype','nucl','-parse_seqids','-out',str(tmp/'vector')],capture_output=True,text=True)
                if proc.returncode: raise RuntimeError(proc.stderr)
                vector_hits=blast_batch(records,tmp/'vector',tmp,args,vector=True) if records else {}
                for hs in vector_hits.values():
                    for h in hs:
                        key=next((token for token in h.chrom.split('|') if token in vector_names),None)
                        if key is None: raise RuntimeError(f'Unexpected vector database identifier: {h.chrom}')
                        h.chrom=vector_names[key]
                vector_meta=[dict(id=n,length=len(s),sha256=hashlib.sha256(s.encode()).hexdigest()) for n,s in vectors]
            all_rows=[];primary=[]
            sequence_by_query=dict(records)
            for item in info:
                if item['line'] in short:
                    row=dict.fromkeys(FIELDS,'');row.update(Species=species,line=item['line'],Chr='-',candidate_id='1',
                        candidate_rank=1,candidate_count=0,mapping_status='input_too_short',junction_status='unresolved',
                        classification='unresolved',Note=f"Input length {item['length']} < min_len {args.min_len}")
                    rows=[row]
                else:
                    qid=item['query_id'];rows=make_rows(item['line'],species,hits.get(qid,[]),vector_hits.get(qid,[]),
                                                      by_chr,args,Path(db).name.lower()=='irgsp1.0',sequence_by_query[qid])
                primary.append(rows[0]);all_rows.extend(rows)
            pairs=check_pairs(args.pairs,primary,args.pair_max_distance) if args.pairs else []
            template=Path(args.template)
            if not template.is_file() and args.template=='./flankingviz.html': template=Path(api.__file__).with_suffix('.html')
            html=template.read_text(encoding='utf-8')
            vendor=Path(__file__).parent/'vendor/jszip.min.js'
            if vendor.is_file():
                html=html.replace('<script src="https://cdn.jsdelivr.net/npm/jszip@3.10.1/dist/jszip.min.js"></script>',
                                  '<script>'+vendor.read_text(encoding='utf-8')+'</script>')
            injection='<script>window.__EMBEDDED_CSV__='+safe_json(csv_text(primary))+';\nwindow.__FV_CANDIDATES__='+safe_json(all_rows)+';</script>'
            if '<body>' not in html: raise ValueError('HTML template has no <body> tag.')
            html=html.replace('<body>','<body>\n'+injection,1)
            metadata=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                engine_sha256=ENGINE_SHA256,
                elapsed_seconds=time.perf_counter()-start,platform=platform.platform(),python=platform.python_version(),
                elapsed_scope='analysis after annotation loading, before export; external benchmark records complete wall time',
                blast_version=subprocess.run(['blastn','-version'],capture_output=True,text=True).stdout.strip(),
                cpu_count=os.cpu_count(),parameters=vars(args),profile=str(profile) if profile else '',db=db,
                reference_provenance=reference_provenance(profile,api),
                inputs=info,vectors=vector_meta,coordinate_convention='1-based adjacent aligned host base; distances are coordinate differences, not intervening-base counts',
                limits='A single flank does not establish insert structure, copy number, vector-backbone absence, or gene function.')
            contents={'.html':html,'.csv':csv_text(primary),'.candidates.csv':csv_text(all_rows),
                      '.json':json.dumps(dict(primary=primary,candidates=all_rows,pairs=pairs),ensure_ascii=False,indent=2),
                      '.run.json':json.dumps(metadata,ensure_ascii=False,indent=2)}
            if pairs: contents['.pairs.csv']=csv_text(pairs,list(pairs[0]))
            elif args.pairs: contents['.pairs.csv']=csv_text([],['pair_id','line_a','line_b','status','distance_bp','screening_distance_bp','interpretation'])
            for suffix,text in contents.items(): (tmp/('result'+suffix)).write_text(text,encoding='utf-8')
            for suffix in contents: os.replace(tmp/('result'+suffix),paths[suffix])
        print(f'[OK] {out}\n     {len(primary)} input records; {len(all_rows)} candidate records.\n     Unresolved junctions remain explicitly labelled.')
    finally:
        lock.unlink(missing_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seqs_dir", default="./seqs", help="Folder containing .seq/.fa/.fasta files")
    ap.add_argument("--profile", default="", help="Species profile YAML (e.g., ./profiles/rice.yaml). If set, it overrides --db/--rap_annot/--rap_pred/--msu.")
    ap.add_argument("--db", default="./Rice/IRGSP1.0", help="BLAST db prefix (makeblastdb -out prefix) [legacy]")
    ap.add_argument("--rap_annot", default="./Rice/rapdb.gff3", help="RAP annotated locus.gff renamed to rapdb.gff3 [legacy]")
    ap.add_argument("--rap_pred", default="./Rice/rapdb_predicted.gff3", help="RAP predicted locus.gff renamed [legacy]")
    ap.add_argument("--msu", default="./Rice/msu_chr_genes.gff3", help="MSU gene-only, chr01.., gene lines only [legacy]")
    ap.add_argument("--template", default="./flankingviz.html", help="Patched FlankingViz HTML template")
    ap.add_argument("--out", default="./out/flankingviz.html", help="Output HTML path")
    ap.add_argument("--min_len", type=int, default=80, help="Minimum sequence length to process")
    ap.add_argument("--vector", default="", help="Optional full-vector FASTA (one/many records) or GenBank")
    ap.add_argument("--junction-side", choices=["auto","start","end"], default="auto", help="Explicit query-side assumption only when no vector boundary is found")
    ap.add_argument("--min-unaligned-tail", type=int, default=20, help="Minimum unaligned query tail for a BLAST-only boundary estimate")
    ap.add_argument("--max-terminal-overhang", type=int, default=10, help="Maximum unaligned bases at the opposite query end")
    ap.add_argument("--min-host-anchor", type=int, default=40)
    ap.add_argument("--min-vector-anchor", type=int, default=20)
    ap.add_argument("--min-identity", type=float, default=80.0)
    ap.add_argument("--min-vector-identity", type=float, default=90.0)
    ap.add_argument("--ambiguity-ratio", type=float, default=0.95, help="Near-best bitscore ratio; heuristic screening threshold")
    ap.add_argument("--max-junction-gap", type=int, default=30)
    ap.add_argument("--max-junction-overlap", type=int, default=12)
    ap.add_argument("--max-hits", type=int, default=100)
    ap.add_argument("--max-hsps", type=int, default=100, help="Maximum HSPs per subject; reaching the limit is reported")
    ap.add_argument("--max-candidates", type=int, default=100, help="Reported loci per input; truncation is explicit")
    ap.add_argument("--batch-size", type=int, default=500)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--pairs", default="", help="Optional CSV with pair_id,line_a,line_b; no pairing inferred from filenames")
    ap.add_argument("--pair-max-distance", type=int, default=10000, help="Screening threshold only; not a structural validation")
    ap.add_argument("--overwrite", action="store_true", help="Explicitly replace files with the requested output prefix")
    args = ap.parse_args()
    for key in ["min_len","min_unaligned_tail","min_host_anchor","min_vector_anchor","max_hits","max_hsps","max_candidates","batch_size","threads"]:
        if getattr(args,key)<1: ap.error(key+" must be positive")
    if not 0 < args.ambiguity_ratio <= 1: ap.error("ambiguity-ratio must be in (0,1]")
    for key in ["min_identity","min_vector_identity"]:
        if not 0 <= getattr(args,key) <= 100: ap.error(key+" must be between 0 and 100")
    for key in ["max_terminal_overhang","max_junction_gap","max_junction_overlap","pair_max_distance"]:
        if getattr(args,key)<0: ap.error(key+" must be nonnegative")

    if args.max_terminal_overhang >= args.min_unaligned_tail:
        ap.error("max-terminal-overhang must be smaller than min-unaligned-tail")

    seqs_dir = Path(args.seqs_dir).resolve()
    if not seqs_dir.exists():
        raise SystemExit(f"[ERROR] seqs_dir not found: {seqs_dir}\nCreate it and put .seq files inside.")

    # Resolve profile / references
    profile_path = Path(args.profile).resolve() if args.profile else None
    if profile_path:
        profile_path = Path(profile_path)

        if profile_path.exists() and profile_path.is_file():
            resolved_profile = profile_path

        else:
            script_dir = Path(__file__).resolve().parent  # = ~/FlankingViz
            name = profile_path.name                      
            candidate = script_dir / "profiles" / f"{name}.yaml"

            if candidate.exists():
                resolved_profile = candidate
            else:
                raise FileNotFoundError(
                    f"Profile '{profile_path}' not found. "
                    f"Expected '{candidate}' or a valid YAML path."
                )

        profile_path = resolved_profile
    gene_maps: List[Dict[str, List[Gene]]] = []

    if profile_path and profile_path.exists():
        cfg = load_profile_yaml(profile_path)

        # Base directory is the profile file's directory (profiles/)
        # Genome/annotation paths can be absolute or relative to the project root.
        # We recommend absolute paths in YAML, but relative paths will be resolved from project root (.. from profiles/).
        project_root = profile_path.parent.parent.resolve()

        def rpath(p: str) -> Path:
            pp = Path(p)
            if pp.is_absolute():
                return pp
            return (project_root / pp).resolve()

        db_setting = str(cfg.get("blast_db", "")).strip()
        db = str(rpath(db_setting))
        if not db_setting:
            raise SystemExit(f"[ERROR] blast_db missing in profile: {profile_path}")

        gff_files = cfg.get("gff_files", [])
        if isinstance(gff_files, str):
            gff_files = [gff_files]
        if not isinstance(gff_files, list) or not gff_files:
            raise SystemExit(f"[ERROR] gff_files missing/empty in profile: {profile_path}")

        # Optional per-profile: feature type for generic parser (default 'gene')
        gene_feature = str(cfg.get("gene_feature", "gene")).strip() or "gene"

        attr_priority = cfg.get("attr_priority", ["Name", "gene", "ID"])
        if isinstance(attr_priority, str):
            attr_priority = [attr_priority]

        for gff_rel in gff_files:
            gff_path = rpath(str(gff_rel))
            if not gff_path.exists():
                raise SystemExit(f"[ERROR] GFF missing: {gff_path} (from {profile_path})")

            bn = gff_path.name.lower()
            if bn == "rapdb.gff3":
                gene_maps.append(load_genes_rap(gff_path, "RAP_ANN"))
            elif bn == "rapdb_predicted.gff3":
                gene_maps.append(load_genes_rap(gff_path, "RAP_PRED"))
            elif "msu" in bn and bn.endswith(".gff3"):
                gene_maps.append(load_genes_msu(gff_path))
            else:
                # Generic GFF3 (Soybean/Maize)
                gene_maps.append(
                    load_genes_generic_gff3(
                        gff_path,
                        source=gff_path.stem,
                        feature=gene_feature,
                        attr_priority=attr_priority,
                    )
                )

        species_name = str(cfg.get("name", "")) or seqs_dir.name
        if not args.vector and cfg.get("vector_fasta"):
            args.vector = str(rpath(str(cfg["vector_fasta"])))

    else:
        # Legacy (rice-only) arguments
        db = str(Path(args.db).resolve())

        rap_annot = Path(args.rap_annot).resolve()
        rap_pred = Path(args.rap_pred).resolve()
        msu = Path(args.msu).resolve()

        if not rap_annot.exists():
            raise SystemExit(f"[ERROR] rap_annot missing: {rap_annot}")
        if not rap_pred.exists():
            raise SystemExit(f"[ERROR] rap_pred missing: {rap_pred}")
        if not msu.exists():
            raise SystemExit(f"[ERROR] msu missing: {msu}")

        gene_maps = [
            load_genes_rap(rap_annot, "RAP_ANN"),
            load_genes_rap(rap_pred, "RAP_PRED"),
            load_genes_msu(msu),
        ]
        species_name = seqs_dir.name

    import sys
    try:
        run_pipeline(args, gene_maps, species_name, db, profile_path, sys.modules[__name__])
    except (ValueError, RuntimeError, OSError) as exc:
        raise SystemExit(f"[ERROR] {exc}")


if __name__ == "__main__":
    main()
