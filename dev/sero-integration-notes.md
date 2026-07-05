# Seroprevalence pathway: integration notes

Status as of this note: seroprevalence support was built and verified against
synthetic data during the jurisdiction-dimension refactor, then deliberately
pulled back out of the MVP. This file exists so that work isn't lost, and so
whoever picks it up next doesn't need to re-derive the design from scratch.
Delete this file once sero is actually re-integrated.

## Why it was pulled

1. No real-data validation exists. Everything sero-related was only ever run
   against fabricated data in a scratch reprex. `tests/test_workflow/testing.R`
   has *always* had a non-functional, commented-out sero block (referencing
   `sero_dat`/`sero_size_mat`/`sero_conversion`, none of which are ever
   assigned) -- there has never been a working real-data example for sero in
   this package.
2. The `discrete_weights` choice for seroconversion (see below) is a real
   modelling decision, not just a software one, and deserves a second pair of
   eyes (e.g. Tianxiao/Nick) before it ships.
3. We didn't want `if`/dispatch logic for a not-yet-validated pathway sitting
   in the MVP's core functions (`fit_waves()`, `stack_jurisdictions()`,
   `prepare_observation_data()`). Rather than merge it in behind conditionals,
   we're waiting until it can go in as a first-class part of the existing
   function workflow.

## The modelling decision: discrete_weights, not discrete_pmf

Case/hospitalisation notification is a one-time event, correctly modelled as
a normalised `discrete_pmf` (delay from infection to that single event).
Seroconversion is different: a person may test positive for many consecutive
days, so the "delay from infection" isn't a probability distribution over a
single event -- it's a persistence/detectability curve that doesn't need to
sum to 1. `epiwave.params::discrete_weights`/`discrete_weights_series` (via
`as_discrete_weights()`) is the right object for this.

`new_convolution_matrix()`/`evaluate()` were widened during this work to
accept `discrete_weights`/`discrete_weights_series` alongside
`discrete_pmf`/`discrete_pmf_series` -- mechanically identical (a
day-difference matrix looked up against a `step` column), just using
`$weight` instead of `$prob`, and *not* forcing normalisation. This was
reverted when sero was pulled (no other current consumer), but re-widening
it is a small, mechanical change -- see the `remove-jurisdiction-dimension`
PR (#51) history for the exact diff if useful as a reference, since the
widening itself was never in question, just its inclusion in the MVP.

## Recommended architecture when re-integrating

Do **not** re-add `define_sero_data()`/`create_small_sero_model()` as
parallel functions with an `if`-based dispatch in `fit_waves()` (that's the
shape we just removed). Instead:

1. Add `total_pop`/`size_vec` as **optional** arguments to
   `define_observation_data()` (default `NULL`), and delete the separate
   `define_sero_data()` wrapper. `prepare_observation_data()` already treats
   `total_pop`/`size_vec` as generically-optional fields (checks
   `'total_pop' %in% names(observation_data)`) -- it doesn't care which
   constructor produced the list, so this needs no changes there once
   `total_pop`/`size_vec` are threaded through.
2. Merge `create_small_sero_model()`'s binomial-likelihood body into
   `create_observation_model()` as an internal branch, dispatched on
   `is.null(observations$total_pop)` -- i.e. move the dispatch that
   currently lived in `fit_waves()` (`model_fn <- if (!is.null(stream$total_pop))
   create_small_sero_model else create_observation_model`) into the body of
   a single `create_observation_model()`, and delete
   `create_small_sero_model()`. `fit_waves()` then calls one function per
   stream uniformly, no dispatch needed there at all.
3. `stack_jurisdictions()`'s `stack_stream()` helper needs its `total_pop`/
   `size_mat` stacking logic re-added (this was mechanical, not
   conceptually hard -- see PR #51 history) -- gated on the same
   `total_pop`-presence check, no new dispatch pattern needed.

This was discussed and agreed as the target shape *before* pulling sero back
out, specifically so this note doesn't need to re-litigate the "one function
vs two" question -- it's already decided, just not implemented yet.

## Real bugs found in the sero pathway during this work (fix again if the old code is reused as a reference)

If anyone copies from the pre-removal sero code (e.g. from PR #51's history)
rather than writing fresh, watch for these -- all were real, silent bugs in
the *original* (pre-refactor) sero code, fixed once during this work:

1. `prepare_observation_data()` stored the seroprevalence sample size as
   `size_vec`, but `create_small_sero_model()` read a field called
   `size_mat` -- always `NULL`. Make sure whatever field name
   `create_observation_model()`'s binomial branch reads matches what
   `stack_jurisdictions()` actually produces.
2. The old `create_small_sero_model()` called
   `new_convolution_matrix(delays, x, n_dates)`, a stale 3-argument call that
   doesn't match the current 2-argument `new_convolution_matrix(pmf, n)`
   signature -- would have errored if ever actually invoked.
3. `fit_waves()`'s call to `create_small_sero_model()` was commented out
   entirely -- the sero pathway was unreachable in the original code, not
   just unvalidated.

## What's still needed before this ships

- A real (or realistic synthetic, committed as an actual test) seroprevalence
  reprex end to end.
- Domain sanity-check on the `discrete_weights` persistence-curve modelling
  choice.
- A committed test, not just a manual script -- there's no automated test
  suite in this package at all yet (`tests/testthat/` isn't populated), so
  this would ideally land alongside seeding that.
