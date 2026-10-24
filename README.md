# Plasmid Verify

A local-first MVP for checking whether a recovered Golden Gate cloning product matches an intended plasmid reference.

## Run

Open `index.html` in a browser. No build step or dependency installation is required.

## Current scope

- Backbone and ordered insert FASTA/raw sequence inputs
- Observed cloned-plasmid FASTA/raw sequence input
- FASTQ/FASTA file upload for the observed sample
- Basic FASTQ read mapping and consensus generation against the expected construct
- Read count and consensus coverage summary
- Expected construct generation
- Golden Gate-style reference matching: backbone + inserts in the intended order
- Optional direct recovered-plasmid versus backbone comparison
- Visual reference-mismatch report with feature map, coordinate ruler, and Sample/Reference sequence window
- Cloning-oriented pass/review result with region-level insert checks
- Full-length plasmid alignment using a memory-efficient banded algorithm (suitable for plasmids around 10 kb)
- Anchor-based alignment for internal insertions, so downstream backbone sequence remains aligned after a large inserted block
- Substitution, insertion, and deletion reporting
- Region-level backbone/insert status
- Circular-plasmid option captured in the interface
- Downloadable text report

Quality-aware FASTQ alignment, quality trimming, coverage plots, circular coordinate rotation, paired-end handling, and long-read assembly are the next production features. The current FASTQ mapper uses a simple approximate match against the expected construct and is a screening aid, not experimental confirmation.

