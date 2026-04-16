"""
Generate ExomeFlow pipeline workflow figure (matches workflow.png style).
Output: workflow_figure.png + workflow_figure.svg
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ── Colours (matched to reference image) ────────────────────────────────────
C_FASTP   = "#F4845F"   # coral/salmon
C_BWA     = "#4A90D9"   # medium blue
C_BAM     = "#8E6BBF"   # medium purple  (SortSam / flagstat / MarkDup)
C_BQSR    = "#6A4C9C"   # deep purple
C_HC      = "#3A7DC9"   # slightly darker blue
C_FILTER  = "#3DAA6E"   # green
C_MERGE   = "#D95B8A"   # pink/magenta
C_SELECT  = "#2A9D8F"   # teal
C_ANNOVAR = "#2C3E7A"   # dark navy
C_ARROW   = "#555555"
C_TEXT    = "white"
C_TITLE   = "#1A1A2E"

FIG_W, FIG_H = 6.5, 18
BOX_W  = 4.6      # main box width
BOX_WS = 2.1      # split box width
BOX_H  = 0.78     # box height
GAP    = 0.22     # gap between boxes
CX     = FIG_W / 2   # centre x

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis("off")
fig.patch.set_facecolor("white")


def box(ax, cx, y, w, h, color, title, lines, radius=0.12):
    """Draw a rounded rectangle with title + bullet lines."""
    x = cx - w / 2
    rect = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=0, facecolor=color, zorder=3,
    )
    ax.add_patch(rect)

    # Title
    n_lines = len(lines)
    title_y = y + h - 0.165
    ax.text(cx, title_y, title,
            ha="center", va="center", fontsize=8.5, fontweight="bold",
            color=C_TEXT, zorder=4)

    # Bullet lines
    if lines:
        line_gap = (h - 0.32) / max(len(lines), 1)
        for i, line in enumerate(lines):
            ly = title_y - 0.19 - i * line_gap
            ax.text(cx, ly, line,
                    ha="center", va="center", fontsize=6.8,
                    color=C_TEXT, zorder=4, linespacing=1.2)


def arrow(ax, cx, y_top, y_bot):
    """Draw a downward arrow."""
    ax.annotate("",
        xy=(cx, y_bot + 0.01),
        xytext=(cx, y_top - 0.01),
        arrowprops=dict(arrowstyle="-|>", color=C_ARROW,
                        lw=1.5, mutation_scale=12),
        zorder=2,
    )


def split_arrow(ax, x_from, y_from, x_to, y_to):
    """Arrow from centre down then across to a split box."""
    ax.annotate("",
        xy=(x_to, y_to),
        xytext=(x_from, y_from),
        arrowprops=dict(arrowstyle="-|>", color=C_ARROW,
                        lw=1.5, mutation_scale=12,
                        connectionstyle="arc,angleA=-90,angleB=90,rad=0"),
        zorder=2,
    )


def merge_arrow(ax, x_from, y_from, x_to, y_to):
    """Arrow from split box back to centre."""
    ax.annotate("",
        xy=(x_to, y_to),
        xytext=(x_from, y_from),
        arrowprops=dict(arrowstyle="-|>", color=C_ARROW,
                        lw=1.5, mutation_scale=12,
                        connectionstyle="arc,angleA=-90,angleB=90,rad=0"),
        zorder=2,
    )


# ── Layout: build top-to-bottom ──────────────────────────────────────────────
y = FIG_H - 0.55

# Raw FASTQ label
ax.text(CX, y, "Raw FASTQ", ha="center", va="center",
        fontsize=11, fontweight="bold", color=C_TITLE, zorder=4)
y -= 0.38
arrow(ax, CX, y + 0.05, y - GAP)
y -= GAP

# Step 1 — fastp
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_FASTP,
    "Step 1: fastp",
    ["Quality Control & Adapter Trimming",
     "• Length ≥ 50 bp,  Base Quality ≥ Q30"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 2 — BWA MEM
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_BWA,
    "Step 2: BWA MEM",
    ["Align to hg38",
     "• -Y  -K 100000000",
     "• Read-group tags set"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 3 — SortSam
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_BAM,
    "Step 3: GATK SortSam",
    ["Coordinate-sort BAM"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 4 — flagstat
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_BAM,
    "Step 4: samtools flagstat",
    ["Alignment QC"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 5 — MarkDuplicates
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_BAM,
    "Step 5: GATK MarkDuplicates",
    ["Remove PCR duplicates"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 6 — BuildBamIndex
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_BAM,
    "Step 6: GATK BuildBamIndex",
    ["BAM Index (.bai)"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 7 — BQSR
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_BQSR,
    "Step 7: GATK BQSR",
    ["Base Recalibration",
     "• dbSNP, Mills, Known Indels",
     "• → recalibrated.bam"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 8 — HaplotypeCaller
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_HC,
    "Step 8: GATK HaplotypeCaller",
    ["Variant Calling",
     "• Exome + dbSNP Annotation"])
y_after_hc = y   # save for split arrows

# ── Variant Filtration label ─────────────────────────────────────────────────
y -= GAP + 0.1
ax.text(CX, y, "Variant Filtration", ha="center", va="center",
        fontsize=8, fontstyle="italic", color=C_ARROW, zorder=4)
y -= 0.18

# Split — Step 9 SNP  |  Step 10 INDEL
y -= BOX_H
cx_snp   = CX - BOX_WS / 2 - 0.12
cx_indel = CX + BOX_WS / 2 + 0.12

# Arrows from HaplotypeCaller down to splits
split_arrow(ax, CX,      y_after_hc - 0.01,  cx_snp,   y + BOX_H + 0.01)
split_arrow(ax, CX,      y_after_hc - 0.01,  cx_indel, y + BOX_H + 0.01)

box(ax, cx_snp,   y, BOX_WS, BOX_H, C_FILTER,
    "Step 9: SNP Filters",
    ["SNP Filtering"])
box(ax, cx_indel, y, BOX_WS, BOX_H, C_FILTER,
    "Step 10: INDEL Filters",
    ["INDEL Filtering"])

y_split_bot = y   # bottom of split boxes

# ── MergeVcfs ────────────────────────────────────────────────────────────────
y -= GAP + 0.18
y_merge_top = y
y -= BOX_H

# Arrows from both split boxes down to MergeVcfs
merge_arrow(ax, cx_snp,   y_split_bot - 0.01, CX, y + BOX_H + 0.01)
merge_arrow(ax, cx_indel, y_split_bot - 0.01, CX, y + BOX_H + 0.01)

box(ax, CX, y, BOX_W, BOX_H, C_MERGE,
    "MergeVcfs",
    ["Merge Variants"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 11 — SelectVariants
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_SELECT,
    "Step 11: SelectVariants",
    ["Extract PASS Variants"])
arrow(ax, CX, y, y - GAP)
y -= GAP

# Step 12 — ANNOVAR
y -= BOX_H
box(ax, CX, y, BOX_W, BOX_H, C_ANNOVAR,
    "Step 12: ANNOVAR",
    ["Functional Annotation",
     "• refGene, ClinVar, gnomAD, COSMIC",
     "• → multianno.vcf + multianno.txt"])

# ── Caption ──────────────────────────────────────────────────────────────────
ax.text(CX, y - 0.32,
        "Whole-exome sequencing variant discovery and annotation pipeline.",
        ha="center", va="center", fontsize=7.5,
        color="#444444", style="italic", zorder=4)

plt.tight_layout(pad=0)
fig.savefig("/media/drprabudh/m2/Exomeflow/workflow_figure.png",
            dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig("/media/drprabudh/m2/Exomeflow/workflow_figure.svg",
            bbox_inches="tight", facecolor="white")
plt.close()
print("Saved: workflow_figure.png  +  workflow_figure.svg")
