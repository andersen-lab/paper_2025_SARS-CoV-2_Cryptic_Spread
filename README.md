# Cryptic Spread of SARS-CoV-2 Undermined Targeted Travel Restrictions

## Delta Variant

- **Time-Scaled Phylogeny:** [`Delta.xml`](./BEAST/Delta/Delta.xml)
- **Discrete Phylogeographic Model:** [`Delta_discrete.xml`](./BEAST/Delta/Delta_discrete.xml)
- **Replicates:** [`replicate0/`](./BEAST/Delta/replicate0/), [`replicate1/`](./BEAST/Delta/replicate1/), [`replicate2/`](./BEAST/Delta/replicate2/)
  - **ESS Log:** `Delta.log`
  - **Jump History:** `jumpHistory.csv.zst` / `jumphistory.parquet`
  - **Persistence Summarizer Log:** `persistence_data/region.states.log`
  - **Persistence Summary:** `persistence_data/summary/`
- **Persistence Summarizer Script:** [`persistence/summarize_persistence_data.py`](./BEAST/Delta/persistence/summarize_persistence_data.py)

## Omicron Variant

- **Time-Scaled Phylogeny:** [`Omicron.xml`](./BEAST/Omicron/Omicron.xml)
- **Discrete Phylogeographic Model:** [`Omicron_discrete.xml`](./BEAST/Omicron/Omicron_discrete.xml)
- **Replicates:** [`replicate0/`](./BEAST/Omicron/replicate0/), [`replicate1/`](./BEAST/Omicron/replicate1/)
  - **ESS Log:** `Omicron.log`
  - **Jump History:** `jumpHistory.csv.zst` / `jumpHistory.parquet`
  - **Persistence Summarizer Log:** `persistence_data/region.states.log`
  - **Persistence Summary:** `persistence_data/summary/`
- **Persistence Summarizer Script:** [`persistence/summarize_persistence_data.py`](./BEAST/Omicron/persistence/summarize_persistence_data.py)

## Supplementary Data

- **Delta Restrictions for investigated countries:** [`delta_restrictions.csv`](./supplementary/delta_restrictions.csv)
- **Omicron Restrictions for investigated countries:** [`omicron_restrictions.csv`](./supplementary/omicron_restrictions.csv)

## Pipeline

Consensus sequences were assembled using this Snakemake pipeline with bwa-mem and iVar v1.2.2.

## Folder Structure

```
.
├── BEAST/
│   ├── Delta/
│   │   ├── Delta.xml
│   │   ├── Delta_discrete.xml
│   │   ├── persistence/
│   │   │   └── summarize_persistence_data.py
│   │   ├── replicate0/
│   │   │   ├── Delta.log
│   │   │   ├── jumpHistory.csv.zst
│   │   │   └── persistence_data/
│   │   │       ├── region.states.log
│   │   │       └── summary/
│   │   ├── replicate1/
│   │   │   ├── Delta.log
│   │   │   ├── jumphistory.parquet
│   │   │   └── persistence_data/
│   │   │       ├── region.states.log
│   │   │       └── summary/
│   │   └── replicate2/
│   │       ├── Delta.log
│   │       ├── jumphistory.parquet
│   │       └── persistence_data/
│   │           ├── region.states.log
│   │           └── summary/
│   └── Omicron/
│       ├── Omicron.xml
│       ├── Omicron_discrete.xml
│       ├── persistence/
│       │   └── summarize_persistence_data.py
│       ├── replicate0/
│       │   ├── Omicron.log
│       │   ├── jumpHistory.csv.zst
│       │   └── persistence_data/
│       │       ├── region.states.log
│       │       └── summary/
│       └── replicate1/
│           ├── Omicron.log
│           ├── jumpHistory.parquet
│           └── persistence_data/
│               ├── region.states.log
│               └── summary/
├── pipeline/
│   ├── scripts/
│   ├── config.json
│   └── Snakefile
└── supplementary/
    ├── delta_restrictions.csv
    └── omicron_restrictions.csv
```

**Data Source:**  
Tegally et al. (2023): [https://github.com/CERI-KRISP/SARS_CoV_2_VOC_dissemination](https://github.com/CERI-KRISP/SARS_CoV_2_VOC_dissemination)
