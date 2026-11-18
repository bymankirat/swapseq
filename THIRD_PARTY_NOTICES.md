# Third-party components

SwapSeq itself is licensed under GPL-3.0-or-later. `index.html` bundles these components, each
under its own licence, all compatible with the GPL.

| Component | Use | Licence |
| --- | --- | --- |
| [JSZip](https://github.com/Stuk/jszip) | Opens QIIME 2 `.qza` artifacts | MIT (dual MIT / GPLv3; used under MIT) |
| [jsfive](https://github.com/usnistgov/jsfive) | Reads BIOM 2.1 (HDF5) tables | Public domain (NIST) |
| [pako](https://github.com/nodeca/pako) | Decompression used by jsfive | MIT and Zlib |
| [Bricolage Grotesque](https://github.com/ateliertriay/bricolage) | Headings | SIL Open Font License 1.1 (`fonts/OFL-bricolagegrotesque.txt`) |
| [Instrument Sans](https://github.com/Instrument/instrument-sans) | Body text | SIL Open Font License 1.1 (`fonts/OFL-instrumentsans.txt`) |
| [JetBrains Mono](https://github.com/JetBrains/JetBrainsMono) | Sample IDs | SIL Open Font License 1.1 (`fonts/OFL-jetbrainsmono.txt`) |

## Data

The demo study inside `index.html` is a subset of the public Student Microbiome Project
(Flores GE et al. 2014, *Genome Biology* 15:531; tables from
[caporaso-lab/student-microbiome-project](https://github.com/caporaso-lab/student-microbiome-project)),
with participant IDs replaced and two labelling errors planted on purpose.

The validation in `validation/` used public tables from the studies listed in `README.md`.
No raw data from those studies is redistributed here.
