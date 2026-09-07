# Validation record

- Profile inventory: 36/36 merged logical conditions present.
- Merged profile rows: 2,100, exactly matching
  `4 boundary conditions x 3 temperature pairs x (25+50+100) sites`.
- Every profile has `samples=512` in its metadata header.
- Non-finite `mean_I` or `se_mean_I` values: 0.
- Regenerated midpoint slopes: 12/12 agree with the frozen analysis to a
  maximum absolute difference of `5.0e-16`.
- Figure outputs: vector PDF plus 600-dpi PNG for both figures.
- LaTeX compile check: PASS under TeX Live/latexmk, with no undefined-reference,
  overfull-box, or underfull-box warning found in the final log.
- Visual PDF check: pages containing the rewritten text, 12-panel profile atlas,
  scaling table, and midpoint figure were rendered and inspected. No clipping,
  overlap, unreadable label, or broken figure was found.

The local compile-check source contains labeled boxes only for the unrelated
`flux_1.png` and `flux2.png` files absent from the downloaded materials. The
revised manuscript source itself retains the original references to those files.

