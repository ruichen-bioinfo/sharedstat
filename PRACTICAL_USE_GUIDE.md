# Practical use: when this package applies, when it only diagnoses, and when it must not be used

This is the guide to reach for before running anything. It answers five questions: whether you have the problem at all,
what your design lets you claim, when the instrumented correction is admissible, which function answers which question,
and when to stop.

This guide describes version 0.1.2, which exports thirteen functions. Version 0.1.0 exports eleven and does not contain
`reliability_from_blocks()` or `reliability_from_replicates()`, so cite the version you actually used, by its own
version DOI rather than by a concept DOI.

## 1. Do you have a shared-component problem?

One test, and it is about the design rather than the statistics:

> Is the same measurement used twice, once inside the quantity you are comparing and once as the reference, baseline or
> predictor you compare it against?

If the answer is no, this package does not apply and forcing it on will produce numbers that mean nothing. Cases where
the answer is yes:

- a contrast of a treated arm against a reference arm, regressed on that same reference arm;
- a ratio whose denominator is also the variable you stratify or order by;
- a residual from a fit whose predictor is measured with the same noise as the response;
- a share of a total, compared across conditions where the total moved.

The error enters twice with opposite signs, so the bias does not shrink as you add features. A larger dataset makes the
estimate more precise and no less wrong.

## 2. What does your design let you claim?

| what you have | what you can do |
|---|---|
| no replication of the reference arm | nothing on the reliability, gain and IV route: reliability cannot be measured, so shared-reference bias cannot be bounded from the data. Four other functions are unaffected -- see below |
| two replicate blocks of the reference | measure reliability, report the naive estimator's expected bias. **Diagnose only** |
| three or more replicates | the same, from every pair, with the spread across pairs as a check |
| an independent second block that shares the signal but not the focal arm's error | the instrumented correction becomes available, conditionally |

The middle two rows are the common case and they are a real result. "Your naive gain of 1.10 is expected to read 0.77 at
the reliability we measured" is a finding, and it needs no correction to be worth reporting.

## 3. When is the instrumented estimator admissible?

Three conditions, all on the design and none testable from the data alone:

1. the instrument shares the latent signal with the reference arm;
2. the instrument does **not** share the focal arm's measurement error: not the same library, not the same replicate,
   not the same batch;
3. the first stage is strong enough to report, and you report it.

`iv_gain()` returns the first-stage strength with the estimate for that reason. If the second condition rests on
assertion rather than on how the samples were prepared, say so and report a sensitivity analysis instead of a corrected
point estimate. A correction whose exogeneity cannot be argued from the design is not a correction; it is a second
estimate with unknown bias.

## 4. Which function answers which question?

| your question | function |
|---|---|
| how reliably was the reference arm measured, from two blocks? | `reliability_from_blocks()` |
| the same, from three or more replicates | `reliability_from_replicates()` |
| where will the naive estimator land, given that reliability? | `naive_gain_bias()` |
| what does the instrumented estimator give, and how strong is its first stage? | `iv_gain()` |
| how large is the signal once the replicate noise floor is removed? | `rms_true()` |
| can I read a count of significant features as an ordering of effect sizes? | `count_vs_amplitude()` |
| is a compositional share non-additive across my contrast? | `share_nonadditivity()` |
| what does a class look like read against a designated reference class? | `reference_class_medians()` |
| cross-block estimate of compensatory divergence, and its attenuation | `crossblock_compensation()`, `attenuation_factor()` |
| what correlation does that imply at minimum? | `compensation_lower_bound()` |
| a known-truth check with compositionality built in | `simulate_closed_sum()` |
| have I met the preconditions at all? | `admission_checklist()` |

`admission_checklist()` is the one to run first. It asks five questions about a planned analysis and answering them
takes less time than discovering the answer afterwards.

## 5. When must you not use this?

Each of these produces a number that looks fine and means nothing:

- **No replication of the reference arm, for the reliability, gain and IV route.** Reliability is unmeasurable, so the
  bias cannot be bounded from the data. Report the design limitation instead. This rules out
  `reliability_from_blocks()`, `reliability_from_replicates()`, `attenuation_factor()`, `naive_gain_bias()`,
  `iv_gain()`, `rms_true()`, `crossblock_compensation()` and `compensation_lower_bound()`.

  It does **not** rule out the rest of the package. `count_vs_amplitude()` needs a noise floor per condition,
  `share_nonadditivity()` needs four genotypes, `reference_class_medians()` needs a class label and a named reference
  class, and `simulate_closed_sum()` needs no data at all. Each states its own admissibility conditions; none of them
  asks for a replicated reference arm. `admission_checklist()` is a reporting helper: it takes whatever of the above you
  were able to compute, and `NA` for the rest, so it works on any design and tells you what your design supports.

  The split is by what a function reads, not by topic. Two of the eight above are on that list for reasons worth
  stating: one needs the same quantity measured in two blocks, and one needs an attenuation factor, which needs a
  reliability. Both are the same requirement in a different guise.
- **A reference arm with no variance.** The reliability is undefined, not 1.
- **An instrument that shares a library, a replicate or a batch with the focal arm.** It shares the error you are trying
  to remove, and the correction will move the estimate in an unknown direction.
- **Ordering conditions by significant-feature counts when their noise floors differ.** `count_vs_amplitude()` exists to
  tell you when that is the situation; when it is, the counts are not comparable however large the difference looks.
- **Class labels that are not mutually exclusive or not defined on one gene universe.** A class whose membership depends
  on two significance thresholds inherits both of them, and its size is then a property of the detection limit.
- **Treating a simulated known truth as validation on real data.** `simulate_closed_sum()` tells you that an estimator
  recovers a truth you constructed. It does not tell you the estimator is right on your data.

## A worked example, end to end

`inst/examples/reproduce_table1.R` is the shortest real use and runs on tables bundled with the package. It measures a
reliability, reads off what a naive estimator is expected to report at that reliability, and rebuilds the exaggeration
factors and the share non-additivity. Run it before running anything on your own data: if its numbers do not come out,
the installation is wrong and nothing else you do will be interpretable.

## What this package will not do for you

It will not tell you whether your reference arm is the right reference. It will not decide whether a corrected estimate
is preferable to a reported bias. It will not make a two-replicate design into a four-replicate one. Those are design
decisions, and the diagnostics exist so that they can be made with a number attached.
