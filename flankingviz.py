#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
flankingviz.py
- Local BLASTN against IRGSP1.0 genome FASTA (your makeblastdb output prefix)
- Nearest gene picking from:
  1) RAP-DB annotated genes (rapdb.gff3; locus.gff)
  2) RAP-DB predicted genes (rapdb_predicted.gff3; locus.gff)
  3) MSU/UGA LOC genes (msu_chr_genes.gff3; gene-only, chr01..)
- Single integrated HTML output (one file) using your FlankingViz interactive viewer:
  - Arrow genes
  - Distance + size
  - F key flip (insertion / genes)
  - Drag x-only, edit labels via double-click

Usage (example):
  python3 flankingviz.py \\
    --seqs_dir ./seqs \\
    --db ./project/IRGSP1.0 \\
    --rap_annot ./project/rapdb.gff3 \\
    --rap_pred ./project/rapdb_predicted.gff3 \\
    --msu ./project/msu_chr_genes.gff3 \\
    --out ./out/flankingviz.html

Notes:
- BLAST options are set to match RAP-DB web defaults you showed:
  - blastn
  - max results 100
  - E-value cutoff 0.1
  - filter low complexity ON (dust)
  - lower case filtering ON (soft masking)
- Insertion coordinate heuristic:
  - take the end of the alignment closer to the query start (border side):
    - if genome strand '+' => insertion_pos = sstart
    - if genome strand '-' => insertion_pos = send
  This mimics flanking reads that start at LB/RB and enter gDNA.
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
    txt = path.read_text(encoding="utf-8", errors="ignore").strip()
    if not txt:
        return ""
    if ">" in txt.splitlines()[0]:
        seq = "".join([l.strip() for l in txt.splitlines() if not l.startswith(">")])
    else:
        seq = "".join([l.strip() for l in txt.splitlines()])
    seq = DNA_RE.sub("", seq).upper()
    return seq

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
    """
    Heuristic for flanking read:
      - query starts at border, so the genomic coordinate closest to query start is the insertion side.
      - If alignment is on '+' strand: border side ~ sstart
      - If '-' strand: border side ~ send (because coordinates decrease)
    """
    return hit.sstart if hit.strand == "+" else hit.send


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
    args = ap.parse_args()

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

        db = str(rpath(str(cfg.get("blast_db", ""))))
        if not db:
            raise SystemExit(f"[ERROR] blast_db missing in profile: {profile_path}")

        gff_files = cfg.get("gff_files", [])
        if isinstance(gff_files, str):
            gff_files = [gff_files]
        if not isinstance(gff_files, list) or not gff_files:
            raise SystemExit(f"[ERROR] gff_files missing/empty in profile: {profile_path}")

        # Optional per-profile: feature type for generic parser (default 'gene')
        gene_feature = str(cfg.get("gene_feature", "gene")).strip() or "gene"

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
                gene_maps.append(load_genes_generic_gff3(gff_path, source=gff_path.stem, feature=gene_feature))

        species_name = str(cfg.get("name", "")) or seqs_dir.name

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

    # Read sequences
    seq_files = sorted([p for p in seqs_dir.iterdir() if p.is_file() and p.suffix.lower() in [".seq",".fa",".fasta",".txt"]])
    if not seq_files:
        raise SystemExit(f"[ERROR] No sequences found in {seqs_dir} (min_len={args.min_len}).\n"
                         f"Put files like: sample1.seq, sample2.seq, ...")

    tmp_dir = Path("./.tmp_flank").resolve()
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Build CSV rows for FlankingViz
    # - line: file stem
    # - species: profile name (if --profile), else seqs_dir name

    rows = []
    header = [
        "Species","line","Chr","insertion_strand",
        "left_id","left_strand","left_dist_bp","left_size_bp",
        "right_id","right_strand","right_dist_bp","right_size_bp",
        "Note"
    ]
    rows.append(",".join(header))

    for fp in seq_files:
        name = fp.stem
        seq = read_any_sequence(fp)
        if len(seq) < args.min_len:
            continue

        hit = blast_best_hit(db=db, seq_name=name, seq=seq, tmp_dir=tmp_dir)
        if hit is None:
            rows.append(",".join([species_name,name,"-","","","","","","","","","","No BLAST hit"]))
            continue

        chr_id = hit.sseqid
        strand = hit.strand
        junction_start = min(hit.sstart, hit.send)
        junction_end = max(hit.sstart, hit.send)
        # Keep a single insertion_pos for display/drag (border side), but distances are computed from junction interval.
        ins_pos = infer_insertion_pos(hit)

        # Gather all genes on that chr from RAP annotated + predicted + MSU
        all_genes = union_genes(gene_maps, chr_id)

        ov, ov_rep = pick_overlap_interval(all_genes, junction_start, junction_end)
        note_bits = [
            f"BLAST: {chr_id}:{junction_start}-{junction_end}({strand}), id={hit.pident:.1f}%, len={hit.length}, e={hit.evalue:.2g}"
        ]
        if ov:
            note_bits.append("Overlap: " + ";".join([f"{g.display_id}({g.source})" for g in ov[:6]]) + ("..." if len(ov)>6 else ""))

        # Determine left/right genes
        left_gene = None
        right_gene = None
        left_dist = None
        right_dist = None

        if ov_rep is not None:
            # mark intragenic as X on left; still compute nearest right outside
            left_gene = ov_rep
            left_dist = "X"
            # pick nearest right outside overlap region (start >= pos)
            nr = nearest_right_interval(all_genes, junction_end)
            if nr:
                right_gene, right_dist = nr[0], nr[1]
        else:
            nl = nearest_left_interval(all_genes, junction_start)
            nr = nearest_right_interval(all_genes, junction_end)
            if nl:
                left_gene, left_dist = nl[0], nl[1]
            if nr:
                right_gene, right_dist = nr[0], nr[1]

        def csv_escape(s: str) -> str:
            s = (s or "").replace('"','""')
            if ("," in s) or ("\n" in s) or ("\t" in s):
                return f"\"{s}\""
            return s

        def gene_fields(g: Optional[Gene], dist) -> Tuple[str,str,str,str]:
            if g is None:
                return ("","","","")
            gid = g.display_id
            gstrand = g.strand
            gsize = str(g.size_bp)
            if dist == "X":
                gdist = "X"
            elif dist is None:
                gdist = ""
            else:
                gdist = str(int(dist))
            return (gid, gstrand, gdist, gsize)

        l_id, l_str, l_dist, l_size = gene_fields(left_gene, left_dist)
        r_id, r_str, r_dist, r_size = gene_fields(right_gene, right_dist)

        rows.append(",".join([
            csv_escape(species_name),
            csv_escape(name),
            csv_escape(chr_id.replace("chr","").lstrip("0") or chr_id),   # show as 1..12 in UI; template doesn't care
            csv_escape(strand),  # insertion strand from BLAST strand (+/-); user can flip in UI
            csv_escape(l_id), csv_escape(l_str), csv_escape(l_dist), csv_escape(l_size),
            csv_escape(r_id), csv_escape(r_str), csv_escape(r_dist), csv_escape(r_size),
            csv_escape(" | ".join(note_bits))
        ]))

    csv_text = "\n".join(rows) + "\n"

    template_path = Path(args.template).resolve()
    if not template_path.exists():
        raise SystemExit(f"[ERROR] template not found: {template_path}\n"
                         f"Copy flankingviz.html next to this script or pass --template.")
    template_html = load_template(template_path)

    html = embed_csv_into_html(template_html, csv_text)

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    csv_out = out_path.with_suffix('.csv')
    csv_out.parent.mkdir(parents=True, exist_ok=True)
    csv_out.write_text(csv_text, encoding='utf-8')
    out_path.write_text(html, encoding="utf-8")

    print(f"[OK] Wrote: {out_path}")
    print(f"     Seqs processed: {len([r for r in rows[1:] if r])}")
    print("     Open the HTML in a browser. (Chrome recommended)")
    print("     Tip: Press 'F' to flip strand for selected gene/insertion.")


if __name__ == "__main__":
    main()
