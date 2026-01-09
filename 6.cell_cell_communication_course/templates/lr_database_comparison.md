# Ligand-Receptor Database Comparison Worksheet

## Purpose
Compare L-R databases to understand coverage differences across CCC tools.

---

## Database Overview

| Database | Tool(s) | Species | Interactions | Last Updated |
|----------|---------|---------|--------------|--------------|
| CellChatDB | CellChat | Human, Mouse | ~2,000 | |
| CellPhoneDB | CellPhoneDB | Human | ~1,500 | |
| OmniPath | LIANA | Multi | ~5,000+ | |
| Consensus | LIANA | Multi | Integrated | |

---

## Interaction Categories

### CellChatDB Categories
- [ ] Secreted Signaling
- [ ] ECM-Receptor
- [ ] Cell-Cell Contact

### CellPhoneDB Categories
- [ ] Secreted
- [ ] Membrane-bound
- [ ] Integrin
- [ ] Receptor-Ligand

---

## Pathway Coverage

| Pathway | CellChatDB | CellPhoneDB | OmniPath |
|---------|------------|-------------|----------|
| TGF-β | | | |
| WNT | | | |
| Notch | | | |
| Chemokine | | | |
| Interleukin | | | |
| Growth factors | | | |

---

## Key Interactions to Check

List specific L-R pairs relevant to your study:

| Ligand | Receptor | In CellChatDB | In CellPhoneDB | In OmniPath |
|--------|----------|---------------|----------------|-------------|
| | | | | |
| | | | | |
| | | | | |

---

## R Code to Explore Databases

### CellChatDB

```r
library(CellChat)
CellChatDB <- CellChatDB.human

# View categories
unique(CellChatDB$interaction$annotation)

# Search for specific gene
CellChatDB$interaction[grep("TGFB1", CellChatDB$interaction$ligand), ]
```

### LIANA Resources

```r
library(liana)

# View available resources
show_resources()

# Get specific resource
lr_pairs <- select_resource("CellChatDB")
```

---

## Notes

[Document any findings about database coverage for your specific analysis]

