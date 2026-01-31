# Adjacency scores in syntenyPlotteR

This branch adds **experimental support for visualising DESCHRAMBLER adjacency scores** alongside Evolutionary History (EH) plots in `syntenyPlotteR`.

The workflow has **two steps**:

1. Convert DESCHRAMBLER output into an adjacency score table (`adjS`)
2. Plot EH blocks + adjacency scores using `draw.eh()`

---
## Cloning the repository (experimental branch)

The adjacency score functionality is currently available on a non-default branch of syntenyPlotteR.

To clone the repository and switch to the experimental branch:

```
git clone https://github.com/Farre-lab/syntenyPlotteR.git
cd syntenyPlotteR
git checkout adjs_score
```


Alternatively, you can clone directly into the branch:

```
git clone -b adjs_score https://github.com/Farre-lab/syntenyPlotteR.git
```

After cloning, you can install and test the package as usual.

---

## 1. Converting DESCHRAMBLER output to `adjS`

### Purpose
DESCHRAMBLER reports adjacency scores between **syntenic fragment (SF) pairs**.  
To visualise these scores, they must be converted into **genomic positions along reconstructed ancestral chromosomes (APCFs)**.

The function `deschrambler_to_adjS()` performs this conversion.

---

### Required DESCHRAMBLER files

| File | Description |
|----|----|
| `Ancestor.APCF` | Reconstructed ancestral chromosome structure (order + orientation of SFs) |
| `Ancestor.ADJS` | Adjacency scores between SF pairs |
| `SFs/block_list.txt` | SF coordinates and IDs (used to compute SF lengths) |

---

### Function

```r
deschrambler_to_adjS(
  apcf_file,
  adjs_file,
  block_list_file,
  out_file = NULL,
  ancestor_name = NULL,
  include_ends = FALSE,
  chr_prefix = ""
)
```

---

### Key arguments

| Argument | Meaning |
|--------|--------|
| `apcf_file` | Path to `Ancestor.APCF` |
| `adjs_file` | Path to `Ancestor.ADJS` |
| `block_list_file` | Path to `block_list.txt` |
| `ancestor_name` | Name of the ancestor (overrides APCF header if provided) |
| `out_file` | Optional output file (tab-separated) |
| `include_ends` | Include 0–SF and SF–0 adjacencies (default `FALSE`) |
| `chr_prefix` | Optional prefix for chromosome IDs |

---

### Example

```r
adjS <- deschrambler_to_adjS(
  apcf_file       = "Ancestor.APCF",
  adjs_file       = "Ancestor.ADJS",
  block_list_file = "SFs/block_list.txt",
  ancestor_name   = "Bovid_ancestor_v3",
  out_file        = "Bovid_ancestor_v3.adjS.txt"
)
```

---

### Output format

The output table (written to file or returned) has **four columns**:

```
ancestor   chr   pos   score
```

- `ancestor`: ancestor name
- `chr`: APCF identifier (e.g. `1`, `2`, `APCF_3`)
- `pos`: genomic position (bp) along the APCF
- `score`: adjacency score (assumed to be in `[0,1]`)

This format is directly compatible with `draw.eh()`.

---

## 2. Plotting EH blocks + adjacency scores

### Purpose
`draw.eh()` now supports plotting adjacency scores in a **separate, narrow panel** next to the EH plot, similar to AGV.

Features:
- Orange → red heatmap for adjacency score
- Grey line connecting scores
- Fixed score scale: **0–1**
- Reference lines at **0.2, 0.5, 0.8**

---

### Required input

| File | Description |
|----|----|
| EH alignment file | Standard `syntenyPlotteR` EH input |
| `adjS` file | Output of `deschramblER_to_adjS()` |

---

### Minimal example

```r
draw.eh(
  output    = "EH",
  chrRange  = c("1", "2"),
  data_file = "my_eh_alignments.txt",
  adj_file  = "Bovid_ancestor_v3.adjS.txt",
  directory = "plots"
)
```

This produces:
- One image per chromosome
- EH blocks faceted by target
- A narrow adjacency score panel on the right

---

### Optional tuning

```r
draw.eh(
  output = "EH",
  chrRange = "1",
  data_file = "my_eh_alignments.txt",
  adj_file = "Bovid_ancestor_v3.adjS.txt",
  adj_panel_width = 0.1,
  directory = "plots"
)
```

---

## Notes & assumptions

- Adjacency scores are assumed to be **between 0 and 1**
- Scores outside this range are automatically clamped
- Adjacency positions correspond to **boundaries between adjacent SFs** along APCFs
- End adjacencies (SF–0) are ignored by default

---

## Recommended testing workflow

1. Run DESCHRAMBLER
2. Convert adjacency scores:
   ```r
   adjS <- deschramblER_to_adjS(...)
   ```
3. Plot:
   ```r
   draw.eh(..., adj_file = "ancestor.adjS.txt")
   ```
4. Compare visually with AGV output

---

## Status

⚠️ These functions are **experimental** and intended for testing and feedback.  
APIs and defaults may change before full integration into `syntenyPlotteR`.

