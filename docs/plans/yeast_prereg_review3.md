# Third, focused review: DRAFT rev 4 of the yeast reciprocity / H1 preregistration

*1 Nov 2026. Reviewer: Claude. Reviewed: `docs/preregistration_reciprocity_yeast.md` (DRAFT rev 4, commit 2293b6c), `analysis/yeast_prep.R`, `analysis/yeastract_convert.R`, `docs/reciprocity_yeast_selections/excluded_pmids.txt`, and the earlier reviews `yeast_prereg_review.md` and `yeast_prereg_review2.md`. Read-only: nothing downloaded, no truth opened, nothing scored, the preregistration not edited. Commands run: `Rscript analysis/yeast_prep.R selftest_gate` in the repo; the same stage on a fixed copy in a scratch directory; a depth-1 clone in a scratch directory; `git merge-base --is-ancestor HEAD HEAD`; the converter's label table on six sample labels. "L" means a line of the rev 4 preregistration, "prep L" a line of `yeast_prep.R` and "conv L" a line of `yeastract_convert.R`.*

**Overall: NO-GO to fix the file as it stands. GO once the PI confirms SF3 and the two blockers below are fixed.** Both blockers are a few lines each and need no data. The SF3 consistency items (S1 below) should be fixed in the same edit. The code items (S2–S4) can follow, but they must be done before step 7.

---

## 1. BL1 (positive fixed-date marker): fixed in the gate, broken in the self-test

- **Gate logic: FIXED.** `FIXED_RE = ^> \*Fixed [0-9]{4}-[0-9]{2}-[0-9]{2}\*$` (prep L51). `prereg_is_fixed()` (prep L76) and `gate_prereg()` (prep L77–81) require the file to be clean against HEAD and to contain the marker. Both truth stages reach this check through `gate_pre_truth()` (prep L377, L397, L464). No line of rev 4 matches the pattern: the header's example is inside backticks and does not start the line (L3), and step 1 (L445–448) contains no literal marker or DRAFT string. Nothing in the text can now lock or unlock the gate.
- **Self-test on the current file: passes.** It printed `current FALSE, fixed_copy TRUE, placeholder_copy FALSE`, and the placeholder `YYYY-MM-DD` is rejected.
- **NEW, BLOCKER (trivial): the mandated self-test fails on the fixed file.** L448 says to run `selftest_gate` on the fixed file, where it "must report 'fixed'". But prep L94–95 stops unless it finds exactly one `> **STATUS: DRAFT` line. Fixing the file removes that line, so the stage stops. I verified this on a fixed copy, which gave `Error: selftest: expected exactly one DRAFT status block line, found 0`. Step 1 therefore cannot be completed as written. **Fix:** when there is no DRAFT line, report `current` and skip the copy tests, or build the copies from whichever header line is present.
- *Notes, no fix required:*
  - The marker is accepted anywhere in the file, even while the DRAFT block is still there.
  - The date is not validated, so `2026-99-99` passes.
  - CRLF line endings or a trailing space make the gate refuse, which is a safe failure.
  - Suggested hardening: require that no `STATUS: DRAFT` header line remains, and that `as.Date()` parses the date.

## 2. BL2 (verdicts and α): mostly fixed, but one ambiguity reopens the 5a/5c boundary

**The verdict table, rebuilt from the text.** Rows are applied in order (precedence).

| Precedence / condition | Verdict | H1 |
|---|---|---|
| Synthetic R check fails (L416, step 5) | none; fix and log | blocked by `gate_pre_truth` (prep L382) |
| Gate F fails or is not computable (L115, L129, L132) | 1 | not run |
| Power gate fails after the fallback, P1 not evaluable, or CYC2008 unobtainable (L195, L210, L408) | 1 | tested if P2/H1 is evaluable (SF3) |
| P2 truth not evaluable (L211) | 1 | "not evaluable" |
| P1 fails the IUT (p ≥ 0.025 or d < 0.005) | 4, plus a regulon-loss line if Holm rejects H_P2-harm | tested |
| P1 passes, P2-NI passes, S1-NI rejected | 2 | tested |
| P1 passes, P2-NI passes, S1 evaluable, S1-NI not rejected | 3a, plus "large-complex cost" if S1-harm p < 0.025 and d ≤ −0.005 | tested |
| P1 passes, P2-NI passes, S1 not evaluable | 3b | tested |
| P1 passes, P2-NI fails, H_P1 not rejected by Holm | 5c | tested |
| P1 passes, P2-NI fails, H_P1 rejected, P2-harm rejected with d ≤ −0.005 | 5b | tested |
| P1 passes, P2-NI fails, H_P1 rejected, P2-harm not rejected | 5a | tested |
| **P1 passes, P2-NI fails, H_P1 rejected, P2-harm p below its Holm level but −0.005 < d₂ < 0** | **unassigned (see below)** | tested |

The rows are mutually exclusive. The claim at L431 is correct: given p₁ < 0.025, if Holm rejects H_P2 at its level, H_P1 then gets 0.025 and is rejected, so "H_P1 not rejected" implies "H_P2 not rejected". The partition is **not** exhaustive, though, and one rule has two readings.

- **NEW, BLOCKER (text, two sentences): the effect-size conditions inside Holm are undefined.**
  - H_P1 is "superiority (with d ≥ +0.005)" (L346), and a loss is "P2-harm rejected under this Holm rule, with d ≤ −0.005" (L349). The text never says whether the d-threshold is part of a member's *rejection*, which drives the step-down, or only of the *claim wording*.
  - The two readings give different verdicts in reachable cases. Example: p_harm = 0.005, d₂ = −0.004 (SE ≈ 0.0015, so P2-NI has z ≈ 0.65 and is not shown), p₁ = 0.020, d₁ = 0.008.
    - If d is part of rejection, H_P2 (the smaller p) is "not rejected", the step-down stops, and the verdict is **5c**.
    - If the step-down uses p-values only, H_P2 is rejected, H_P1 is rejected at 0.025, and the verdict is 5a or nothing. "P2-harm is not [rejected]" in 5a (L422) is false, and 5b (L423) fails on d, which gives the unassigned row above.
  - The same ambiguity decides the verdict-4 loss line. Example: p₁ = 0.001 with d₁ = 0.004, and p_harm = 0.02 with d₂ = −0.006.
  - **Fix:**
    - "Holm's step-down uses p-values only. The d-thresholds qualify only the wording of the claim."
    - 5a's condition becomes "and no significant regulon loss (as defined at L349)".
    - Use one inequality convention (p < level) throughout.
- **α accounting: coherent.**
  - Adoption is an IUT: each component at 0.025, so the adoption error is ≤ 0.025.
  - The S family is a serial gatekeeper behind the IUT: S1 at 0.025, then Holm over S2–S6 at 0.025 only if S1-NI is rejected. That is valid.
  - The NI and harm nulls are disjoint (L344, L350), so S1-harm and P2-harm at full level add no FWER.
  - The composite member H_P2 (rejected if NI or harm is rejected at the member's level, with p ordered by the minimum) is a valid local test, because at most one of its nulls is true.
  - Family H is separate and in effect two-sided at 0.05, as disclosed (L281).
  - *Note:* claims made via "IUT if both pass, else Holm" control FWER over {H_P1, H_P2} at α under independence (0.0250 exactly) or positive dependence. Under adversarial negative dependence the worst case is 1.5α. One sentence of disclosure is enough; no change is needed.
- *Note:* the precedence (no verdict, then 1, then the rest) is implied but never stated. Add one line saying so.

## 3. SF2 (git-ancestry enforcement): the idea is right, but the code can be bypassed

`gate_after(later, earlier)` (prep L387–392) checks `merge-base --is-ancestor last_commit(earlier_f) last_commit(later)`.

- **(a) Same-commit bypass (SHOULD-FIX, high).** Git treats a commit as its own ancestor; `git merge-base --is-ancestor HEAD HEAD` exited 0. So committing the 8a counts and the export log *together* passes. That loses the very evidence the rule exists for: that 8a was committed before the export. **Fix:** require `ce != cl` (strict ancestry).
- **(b) The code uses the last commit, not the adding commit (SHOULD-FIX, high).** L148 says "the commit that **adds** the export log", but the code uses the *last* commit touching it (prep L74, L388). A log committed early, even before gate F, and touched again after 8a passes. **Fix:** for `later`, use the first commit that added the file (`git log --diff-filter=A --format=%H -- f | tail -1`), and also require that the log's last commit equals its first, or log the difference as a deviation. Using the last commit for the `earlier` files is correct.
- **(c) A shallow clone passes everything (SHOULD-FIX).** I checked this on a depth-1 clone. Every path's `git log -1` returns HEAD (2293b6c), including `yeast_prereg_review.md`, whose real last change is fd9c827. Every ancestry check then compares HEAD with HEAD and passes. **Fix:** stop if `git rev-parse --is-shallow-repository` returns `true`.
- **(d) Amended or rebased history.** Rewriting history can fabricate the order; the gate fails only if commits are missing. Ancestry also proves the order of commits, not when the YEASTRACT web export was performed. L148 ("performed only after") cannot be enforced by code; say so. **Recommended:**
  - push each gate commit (4, 6, 7, 8a) before the next step;
  - have the gate require each `earlier` commit to be reachable from the remote-tracking branch, so the remote's push record gives an external timestamp;
  - record the raw export file's SHA-256 in the export log.
- **(e) The converter-fix path conflicts with the gate.** L160 allows a logged converter fix "before any association row is used". Such a fix necessarily comes after the export, so the converter's last commit follows the export log's commit, and `gate_after` (prep L466) refuses. The only way through is to re-commit the log, which is loophole (b). **Fix:** say how a converter deviation passes the gate, for example a committed `converter_deviation.md` whose commit must follow the log, while the 8a-recorded `converter_blob` stays the reference.
- **(f) The text says 8c runs the converter; the code does not (SHOULD-FIX, high).** L464 says it "runs the committed converter". `regulon_structure` instead accepts any normalised TSV (`--yeastract`, prep L467–470), with no link to the logged raw export. A hand-edited TSV would pass. **Fix:** 8c should run `yeastract_convert.R` itself on the raw export file(s) and check their SHA-256 against the committed export log.
- **(g) `analysis/yeast_prep.R` is itself not gated (SHOULD-FIX, high).** `need_clean` (prep L378–380) covers the test script and the converter, but not the script that builds both truths. An uncommitted, edited copy can run 8a and 8c, and no blob of it is recorded. **Fix:** add `analysis/yeast_prep.R` to `need_clean` and record its blob in both count files.
- *Note (CYC2008, prep L398, L436–442):*
  - The log is checked by ancestry only when `--cyc2008` is passed. It is not checked with `need_clean`, and the file's hash is not compared with the log.
  - The two-run procedure is not written down. The first run without the file gives "CYC2008 required but not available"; then the file is downloaded and logged, and 8a is re-run, which overwrites the counts.
  - Say that the first run's counts are committed too.

## 4. SF3 (provisional) consistency across the file

These places agree: L27, L48, L195, L210, L279, L406, L429 and L462. Each says that the power gate comes from complex sizes before any regulon membership (the code agrees: `complex_structure`, prep L396–460, never reads YEASTRACT), that the reciprocity arms get verdict 1, and that H1 is still tested if gate F passed and P2/H1 is evaluable, which consumes the regulon truth. **Contradictions to fix together with the PI's confirmation (SHOULD-FIX, high):**

1. **Is R scored?** L408 says "no reciprocity score is computed", but the same paragraph says "the R vs K10 regulon contrast is then reported descriptively". Step 10 (L466) also says "H1 only, with the R vs K10 regulon contrast descriptive". Scoring R against yeast regulons is a reciprocity score on a yeast truth. Pick one. I recommend no R score at all under verdict 1, which matches "synthetic-only candidate" (L417) and "declared the last" (L39).
2. **Human and descriptive arms under verdict 1 from step 8a are unspecified.** L466 ("…then HepG2 and Jurkat, then the descriptive sets") reads as if they run. L408 and L417 suggest they do not.
3. **What "consumed" covers.** L46 says "once the truth phase runs, yeast is consumed for good", which is the whole compendium. L48 says H1 "consumes the yeast regulon truth; the complexes stay unscored", which reads per truth. State whether the unscored complexes stay available for any future test. I recommend not: the deleteome is consumed as soon as any yeast truth is scored.
4. **The gene-jackknife RNG stream.** One seed per dataset, with truths drawn in the order P1, P2, S1, S6 (L330). Whether P1's halves are still drawn under verdict 1 changes P2's and H1's halves. **Fix:** "P1's halves are drawn (and discarded) in every case."
5. *Minor:* the verdict 1 row (L417) says "Otherwise H1 is tested" without "if the regulon truth is evaluable" (compare L279 and L429).

## 5. Other new findings in rev 4

- **Converter label table misclassifies (SHOULD-FIX, high; code; must be fixed before step 7).** In conv L27 the binding regex includes `direct`, which matches **"Indirect"**. Because the entries are checked in order, "Indirect evidence" maps to **binding** (verified). "Expression; ChIP-chip" also maps to binding only, because the combined-label regex needs "and|plus|&". If YEASTRACT uses Direct/Indirect labels, expression-only associations (possibly from the deleteome) would enter the binding leg of the primary truth, and the binding-only fallback would be contaminated too. That is the B1 failure through the back door. **Fix:**
  - use word boundaries, and test "indirect" before "direct" (or use an exact label table);
  - stop if a label matches more than one class unintentionally;
  - **stop on any unmapped label.** At present such rows are silently dropped (conv L39) and only logged (conv L44–45).
- *Note:* the value-column regex `(regulation|type|value|sign)` (conv L38) can pick up an "Evidence type" column. That is harmless for associations, because any non-zero value counts, but the logged value table would then be meaningless.
- *Note:* `excluded_pmids.txt` still says `truth_structure` for the screen additions; the stage is now `regulon_structure`.
- *Note:* nothing checks that `excluded_pmids.txt` was committed at or before the commit that fixes the file, as L143 requires. A strict ancestry check against the preregistration's fixed commit would close this.

---

## Summary

**BLOCKERS (text or one-line code; before fixing the file)**

1. BL2 residual: the d-thresholds inside the Holm step-down are undefined, which leaves a reachable cell unassigned (5a vs 5b) and lets two readings flip 5a and 5c. Fix: step-down on p only, d qualifies the wording, and 5a means "no significant regulon loss (L349)".
2. `selftest_gate` stops on the fixed file (prep L94–95), so step 1 (L448) cannot be completed.

**SHOULD-FIX (high)**

- SF3 consistency: whether R is scored under verdict 1 (L408 vs L466), whether the human and descriptive arms run, what "consumed" covers, and the RNG stream.
- SF2 code:
  - strict ancestry;
  - the adding commit rather than the last commit for the log;
  - refuse shallow clones;
  - put `yeast_prep.R` under `need_clean`;
  - 8c runs the converter on the hashed raw export;
  - a defined route for a converter deviation.
- Converter "Indirect" → binding bug; stop on unmapped labels.

**NOTES:** marker and date hardening; IUT + Holm FWER worst case under negative dependence; precedence line; CYC2008 log and re-run; stale stage name; ancestry of `excluded_pmids.txt`; push gate commits to the remote for external timestamps.

**Go/no-go:** NO-GO as is. **GO to fix the file once the PI confirms SF3 and the two blockers and the SF3 consistency items are edited in.** All of it is about 10 lines of text plus a one-line self-test fix. The SF2 and converter code items need no change to the preregistration text, since they bring the code in line with L148 and L464, and they can follow the fix. They must pass before step 7 commits the converter and the test script. Nothing found affects work done so far: no truth has been opened.
