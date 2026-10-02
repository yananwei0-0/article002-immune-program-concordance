# Article 002 release audit

Technical validation: **PASS on 2026-09-29**
Public release ready: **NO — anonymity, named authors, DOI, and final rerun remain**

Validated on 2026-09-29 against the locked consolidated workbook.

- 23 workbook sheets: one contents sheet and 22 analytical sheets.
- 22 CSV exports matched the workbook cell-for-cell and passed SHA-256 checks.
- Four main and six supplementary figures were present in PNG, SVG, and TIFF formats.
- PNG files passed the 300-dpi check; RGB TIFF files passed the 600-dpi check.
- The regenerated manifest used portable relative paths and the exact locked package versions.
- The manuscript claim was corrected to 77/77 positive RNA–protein correlations and 76/77 with q < 0.05.
- Unsupported acetylation and prespecification wording was removed.
- No absolute local paths or credential-like strings were detected in public text assets.

Metadata preparation update on 2026-10-02:

- Added the scoped MIT/CC BY 4.0 dual-license files and explicit exclusions.
- Confirmed that none of the 105 upstream source assets is present; every source-manifest row is marked `NOT_REDISTRIBUTED`.
- Added the final repository URL and a source-family terms audit.
- Replaced the citation template with a CFF 1.2.0 schema-valid interim `CITATION.cff` using a collective author label.
- This metadata-only update has not yet been followed by a full repository validator rerun.

The remaining blockers are listed in `PUBLIC_RELEASE_BLOCKERS.md` and do not require scientific recomputation.
