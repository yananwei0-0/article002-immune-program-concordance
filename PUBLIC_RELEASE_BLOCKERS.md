# Public-release blockers

The technical reproducibility checks pass. Metadata that does not compromise
review anonymity has been prepared; the following checklist controls the final
public release:

- [ ] Replace author, affiliation, corresponding-author, funding, competing-interest, contribution, and acknowledgement placeholders in the manuscript.
- [x] Approve and add the scoped dual license: MIT for code and CC BY 4.0 for original figures/tables/documentation, with the exclusions in `LICENSE.md`.
- [x] Confirm the package redistribution boundary: all 105 upstream source assets are absent and marked `NOT_REDISTRIBUTED`; source-family terms are audited in `provenance/SOURCE_TERMS_AUDIT.md`.
- [x] Add the final GitHub repository URL to repository and citation metadata.
- [ ] Insert the public repository URL into the manuscript only when journal anonymity rules permit.
- [ ] Create the immutable archived release and add its DOI.
- [x] Replace `CITATION.cff.template` with a schema-oriented interim `CITATION.cff` using the collective author label.
- [ ] Replace the collective citation author with the final named author list and add the archived DOI before public release.
- [x] Run `python validate_release.py` from the complete repository after the metadata changes and review the regenerated checksums.
- [ ] Set `public_release_ready` to `true` only after all preceding items are complete.

No scientific calculation or figure-generation blocker remains in this candidate.
