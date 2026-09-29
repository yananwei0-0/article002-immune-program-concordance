# Public-release blockers

The technical reproducibility checks pass, but the following author-controlled items must be completed before publication:

- [ ] Replace author, affiliation, corresponding-author, funding, competing-interest, contribution, and acknowledgement placeholders in the manuscript.
- [ ] Select and approve the code/documentation/figure license; replace `LICENSE_REQUIRED.md` with the final `LICENSE` file.
- [ ] Confirm redistribution compatibility for every source-manifest row still marked `NOT_ASSESSED`.
- [ ] Add the final GitHub repository URL to the manuscript and citation metadata.
- [ ] Create the immutable archived release and add its DOI.
- [ ] Replace `CITATION.cff.template` with a complete `CITATION.cff`.
- [ ] Run `python validate_release.py` after the metadata changes and review the regenerated checksums.
- [ ] Set `public_release_ready` to `true` only after all preceding items are complete.

No scientific calculation or figure-generation blocker remains in this candidate.
