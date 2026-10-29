#!/usr/bin/env python3
"""Combine every gene_presence_table_v2.csv into one workbook, a tab per genome.

  python3 build_workbook.py <results-dir> <out.xlsx>
"""

import csv
import glob
import os
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo

csv.field_size_limit(2**31 - 1)          # MULTISPECIES nr deflines

FONT = "Arial"
HDR_FILL = PatternFill("solid", fgColor="1F3864")
HDR_FONT = Font(name=FONT, size=10, bold=True, color="FFFFFF")
BODY = Font(name=FONT, size=10)

# final_call -> fill. Only the call column is coloured; the point is to make
# 'not present' visually distinct from 'present' without implying that the
# four not-present states mean the same thing.
CALL_FILL = {
    "present":            PatternFill("solid", fgColor="C6EFCE"),  # green
    "present_frameshift": PatternFill("solid", fgColor="FFEB9C"),  # amber
    "flagged":            PatternFill("solid", fgColor="FFD9CC"),  # orange
    "unverifiable_short": PatternFill("solid", fgColor="E4DFEC"),  # grey-violet
    "absent":             PatternFill("solid", fgColor="F2F2F2"),  # grey
}
CALL_FONT = {
    "present":            Font(name=FONT, size=10, bold=True, color="006100"),
    "present_frameshift": Font(name=FONT, size=10, color="9C6500"),
    "flagged":            Font(name=FONT, size=10, color="9C3314"),
    "unverifiable_short": Font(name=FONT, size=10, color="403151"),
    "absent":             Font(name=FONT, size=10, color="808080"),
}

# width, wrap
COLW = {
    "gene": 10, "round1_call": 15, "round1_identity_pct": 11,
    "round1_coverage_pct": 11, "contig": 20, "coordinates": 17,
    "nr_confirming_accession": 21, "nr_confirming_title": 46,
    "nr_confirming_organism": 30, "nr_identity_pct": 10,
    "nr_hit_coverage_pct": 11, "nr_coverage_basis": 12, "nr_hit_rank": 9,
    "frame_status": 19, "frameshift_at": 13, "nr_verification": 60,
    "final_call": 19,
}
NUMERIC = {"round1_identity_pct", "round1_coverage_pct",
           "nr_identity_pct", "nr_hit_coverage_pct", "nr_hit_rank"}


def genome_dirs(root):
    out = []
    for p in sorted(glob.glob(os.path.join(root, "*", "gene_presence_table_v2.csv"))):
        out.append((os.path.basename(os.path.dirname(p)), p))
    return out


def gene_order(root, samples):
    """Gene order from the wide matrix header, so tabs read like the matrix."""
    wide = os.path.join(root, "summary", "gene_presence_matrix_wide.tsv")
    if os.path.exists(wide):
        with open(wide) as fh:
            return next(csv.reader(fh, delimiter="\t"))[1:]
    with open(samples[0][1]) as fh:
        return [r["gene"] for r in csv.DictReader(fh)]


def write_sheet(ws, header, rows):
    ws.append(header)
    for c in range(1, len(header) + 1):
        cell = ws.cell(row=1, column=c)
        cell.fill, cell.font = HDR_FILL, HDR_FONT
        cell.alignment = Alignment(horizontal="center", vertical="center",
                                   wrap_text=True)
    ws.row_dimensions[1].height = 30

    call_i = header.index("final_call") + 1
    for r in rows:
        ws.append([r.get(h, "") for h in header])
        i = ws.max_row
        for c, h in enumerate(header, start=1):
            cell = ws.cell(row=i, column=c)
            cell.font = BODY
            cell.alignment = Alignment(vertical="top",
                                       wrap_text=h == "nr_verification")
            if h in NUMERIC and cell.value not in ("", None):
                try:
                    cell.value = float(cell.value)
                    cell.number_format = "0.0" if h.endswith("_pct") else "0"
                except ValueError:
                    pass
        call = r.get("final_call", "")
        cc = ws.cell(row=i, column=call_i)
        cc.fill = CALL_FILL.get(call, PatternFill())
        cc.font = CALL_FONT.get(call, BODY)
        cc.alignment = Alignment(vertical="top", horizontal="center")

    for c, h in enumerate(header, start=1):
        ws.column_dimensions[get_column_letter(c)].width = COLW.get(h, 14)
    ws.freeze_panes = "B2"
    ref = f"A1:{get_column_letter(len(header))}{ws.max_row}"
    t = Table(displayName=f"t_{ws.title.replace('-', '_')}", ref=ref)
    t.tableStyleInfo = TableStyleInfo(name="TableStyleLight1",
                                      showRowStripes=False)
    ws.add_table(t)


def main():
    root, out = sys.argv[1], sys.argv[2]
    samples = genome_dirs(root)
    if not samples:
        sys.exit(f"no gene_presence_table_v2.csv under {root}")
    genes = gene_order(root, samples)

    wb = Workbook()
    wb.remove(wb.active)

    # ── overview tab ────────────────────────────────────────────────────────
    # 29 tabs is too many to navigate blind, so the first sheet is the wide
    # matrix: genomes down, genes across, 1/0. The per-genome tabs carry the
    # evidence behind each cell.
    ov = wb.create_sheet("Overview")
    tables = {}
    for s, p in samples:
        with open(p) as fh:
            tables[s] = {r["gene"]: r for r in csv.DictReader(fh)}

    ov.append(["genome"] + genes + ["n_present"])
    for s, _ in samples:
        row = [s]
        for g in genes:
            r = tables[s].get(g)
            row.append("" if r is None else
                       (1 if r["final_call"] == "present" else 0))
        row.append(None)                     # filled with a formula below
        ov.append(row)
    last_gene_col = get_column_letter(1 + len(genes))
    total_col = get_column_letter(2 + len(genes))
    for i in range(2, ov.max_row + 1):
        ov[f"{total_col}{i}"] = f"=SUM(B{i}:{last_gene_col}{i})"
    # A genuinely blank row between the data and the totals. ws.append([])
    # does not advance max_row, so the totals landed directly under the last
    # genome -- one row from being copied into an analysis as a 30th sample.
    tot = len(samples) + 3
    ov.cell(row=tot, column=1, value="n_genomes")
    for c in range(2, 2 + len(genes)):
        L = get_column_letter(c)
        ov.cell(row=tot, column=c, value=f"=SUM({L}2:{L}{len(samples) + 1})")
    ov.cell(row=tot, column=2 + len(genes),
            value=f"=SUM(B{tot}:{last_gene_col}{tot})")

    for c in range(1, 3 + len(genes)):
        cell = ov.cell(row=1, column=c)
        cell.fill, cell.font = HDR_FILL, HDR_FONT
        cell.alignment = Alignment(horizontal="center", vertical="bottom",
                                   text_rotation=90)
        ov.cell(row=tot, column=c).font = Font(name=FONT, size=10, bold=True)
    ov.row_dimensions[1].height = 92
    ov.column_dimensions["A"].width = 13
    for c in range(2, 2 + len(genes)):
        ov.column_dimensions[get_column_letter(c)].width = 4.6
    ov.column_dimensions[total_col].width = 11
    thin = Side(style="thin", color="BFBFBF")
    for row in ov.iter_rows(min_row=2, max_row=len(samples) + 1,
                            max_col=2 + len(genes)):
        for cell in row:
            cell.font = BODY
            cell.alignment = Alignment(horizontal="center")
            cell.border = Border(left=thin, right=thin, top=thin, bottom=thin)
            if cell.value == 1:
                cell.fill = CALL_FILL["present"]
            elif cell.value == 0:
                cell.fill = CALL_FILL["absent"]
        row[0].alignment = Alignment(horizontal="left")
    ov.freeze_panes = "B2"

    # ── legend ──────────────────────────────────────────────────────────────
    # 1/0 is not the whole story: four distinct call states all score 0, and
    # 'flagged' in particular is unresolved rather than absent. Anyone reading
    # the matrix without that will read 0 as "the gene is not there".
    lg = tot + 2
    ov.cell(row=lg, column=1, value="Legend").font = Font(
        name=FONT, size=10, bold=True)
    notes = [
        ("1", "present — nr-confirmed", "present"),
        ("0", "present_frameshift — gene found, intact protein not supported",
         "present_frameshift"),
        ("0", "flagged — locus found, could not be verified. NOT evidence of absence",
         "flagged"),
        ("0", "unverifiable_short — too short for nr to verify either way",
         "unverifiable_short"),
        ("0", "absent — no locus, or below the superfamily threshold", "absent"),
    ]
    for i, (score, text, key) in enumerate(notes, start=1):
        c = ov.cell(row=lg + i, column=1, value=score)
        c.fill, c.font = CALL_FILL[key], CALL_FONT[key]
        c.alignment = Alignment(horizontal="center")
        t = ov.cell(row=lg + i, column=2, value=text)
        t.font = BODY
        t.alignment = Alignment(horizontal="left")
    ov.cell(row=lg + len(notes) + 2, column=1,
            value="Per-genome tabs carry the evidence behind every cell. "
                  "See CAVEATS.md before interpreting sesI, bhp, aap or the "
                  "presence-is-not-resistance genes (norA, dfrC, uhpT).").font = \
        Font(name=FONT, size=10, italic=True)

    # ── one tab per genome ──────────────────────────────────────────────────
    with open(samples[0][1]) as fh:
        header = next(csv.reader(fh))
    for s, p in samples:
        with open(p) as fh:
            rows = list(csv.DictReader(fh))
        rows.sort(key=lambda r: genes.index(r["gene"])
                  if r["gene"] in genes else len(genes))
        write_sheet(wb.create_sheet(s[:31]), header, rows)

    wb.save(out)
    print(f"wrote {out}: {len(wb.sheetnames)} sheets "
          f"({len(samples)} genomes x {len(genes)} genes)")


if __name__ == "__main__":
    main()
