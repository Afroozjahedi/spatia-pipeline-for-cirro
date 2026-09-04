# `spatia/analysis/channels.py`

**Added 2026-09-04.** Not a pipeline step of its own — a support module `segmentation.py` (and, going forward, `preprocessing.py`) calls into for channel-structure detection, panel resolution, and stack normalization.

**Why this exists:** see `segmentation.py`'s 2026-09-04 changelog for the full incident. Short version — the real CRC TMA cores are CODEX/PhenoCycler raw output: ImageJ hyperstacks with axes `TCYX`, shape `(23, 4, 1440, 1920)` (23 acquisition cycles × 4 channels/cycle = 92 planes), not a flat 92-channel stack. That single fact broke two things on 2026-09-03: the embedded per-image labels are scanner placeholders (`"ch1".."ch4"`, no marker identity), so a name-for-name comparison against the real panel could never pass; and `spacec` reads the file itself and always treats axis 0 as the channel axis, so it saw "23 channels" and crashed with an `IndexError` on a perfectly good file.

## What can and cannot be automated

Automatable, because the file's own structure tells us:
- axis layout and plane count
- cycle count and channels-per-cycle
- which planes are the nuclear stain (channel 1 of every cycle, CODEX convention)
- which plane is the last marker (the last one)
- a positional template panel (`"cyc01_ch1"`, ...) when no panel file is supplied

**Not automatable:** the biological identity of each marker. That's experiment metadata from the acquisition sheet, not something present in the pixels — a pipeline that guessed it would be inventing data. The panel file (`channelnames.txt`) stays the one required manual input; this module's job is validating it against the real file (by length, and by a pixel-level periodicity check) rather than trusting it blindly.

## Entry points

```python
from spatia.analysis import channels as ch

info   = ch.inspect_stack(tiff_path)                 # read structure only, no pixels
panel  = ch.ensure_panel(channel_file_path, info)     # validate/create the marker panel
stack  = ch.load_flat_stack(tiff_path, info)          # (C, Y, X) array, acquisition order
resolved = ch.resolve_channels(panel, info,
               nuclei_channel="auto",
               membrane_channel_list=["CD45"])
qc = ch.verify_nuclear_periodicity(stack, info)       # pixel-level order check
ch.normalize_stack_to_file(tiff_path, out_path, panel, info)
```

## Key functions

| Function | Role |
|---|---|
| `inspect_stack(path)` → `StackInfo` | Reads only the TIFF header (cheap — safe to call per-image). Determines axis layout (`YX`, `CYX`, `TCYX`, `ZCYX`), true plane count, cycle structure, and where channel names come from (ImageJ `Labels`, OME-XML, or QPTIFF `ImageDescription`). Raises on any other axis order rather than guessing — a `CTYX` file has the same shape family as `TCYX` but the opposite memory order, and reshaping it as though cycle-major would silently permute every marker. |
| `load_flat_stack(path, info=None)` | Loads pixels and flattens to `(C, Y, X)` in acquisition (cycle-major) order: cycle 1's 4 channels, then cycle 2's, etc. |
| `read_panel(panel_path)` | One channel name per line, in stack order. |
| `panel_template(info)` | Positional placeholder names (`"cyc01_ch1"`, ...) derived purely from structure — honest names that say *where* a plane sits without pretending to know what it stains. |
| `ensure_panel(panel_path, info, auto_template=True)` | Returns the panel for this stack. **The length check here is the validation with real teeth** — it catches a truncated acquisition, a panel from a different experiment, or a core genuinely missing cycles, which is exactly what the old name-equality check was reaching for and missing. Writes a template (with a loud warning) if the file is missing and `auto_template=True`. |
| `classify_panel(names, info)` | Splits panel positions into `nuclear` / `blank` / `marker` index groups. `blank` indices matter downstream: those planes are identically zero, so any per-channel standard-deviation normalization divides by zero there. |
| `resolve_channel(spec, names, info, role)` | Turns one config value into an exact panel name. Accepts `"auto"` (nuclei → channel 1 of cycle 1; last_marker → final plane), a short name matched uniquely by exact/prefix/substring against the panel (`"CD45"` → `"CD45 - hematopoietic cells"`), an exact string, or a zero-based index. **Ambiguity raises** — e.g. `"CD4"` prefix-matches CD44/CD45/CD45RA/CD45RO in the CRC panel; picking the first silently would be a coin flip nobody would notice was wrong. |
| `resolve_channels(names, info, nuclei_channel, membrane_channel_list, last_marker)` | Resolves all three config values in one call, one error style. |
| `verify_nuclear_periodicity(arr, info, min_corr=0.5)` | Pixel-level QC for the flattening order. In CODEX, channel 1 of every cycle re-images the same nuclei — so after a correct flatten, those planes must correlate with each other. If they don't, the flatten order or the panel doesn't match this acquisition. Automates the manual correlation check run by hand on `TMA_A/reg011_X01_Y01_Z09.tif` during this session (result: r = 0.87–0.92 across all 23 nuclear positions + DRAQ5). |
| `normalize_stack_to_file(src_path, dst_path, names, info=None, overwrite=False)` | Writes a flat `(C, Y, X)` ImageJ TIFF with the real marker names embedded in `Labels`, so any downstream tool — `spacec` included — reads the stack the way the panel describes it. Written as plain `.tif`, deliberately **not** `.ome.tif` — `segmentation.py`'s file glob excludes `.ome.tif`, so an OME export here would be silently skipped and the run would report zero files. |

## OME-TIFF / QPTIFF support (added same day, before any real file of either format was run through this pipeline)

`inspect_stack()` also reads channel names from OME-XML (`<Channel Name="...">`, parsed with the stdlib XML parser rather than a regex, since `<Image Name="...">` and `<Plate Name="...">` also carry `Name` attributes a regex would happily collect too) and from QPTIFF per-page `ImageDescription` XML (`<Biomarker>` / `<Name>`).

**The distinction that matters:** `StackInfo.labels_unverifiable` is `True` when a format that *normally* carries channel names (OME-TIFF, QPTIFF) yields none on parsing — this is different from `labels_are_generic()`, which means the file legitimately has no marker identity (raw CODEX). `labels_are_generic()` is explicitly `False` whenever `labels_unverifiable` is `True`. Collapsing these two states would mean applying the panel positionally to a file that may already have real, different channel names embedded — silently mislabeling every marker. `segmentation.py`'s channel check fails loud on the unverifiable case rather than falling through to the count-only comparison used for genuinely generic files.

**Verification status (confidence: high for OME, low for QPTIFF).** OME-TIFF name reading, unnamed-OME detection, and pyramidal-OME handling (3 levels) were unit-tested against synthetic files in this session. **No real `.qptiff` file was available to test against** — the QPTIFF reader is untested against real Akoya/PhenoImager output and should be treated as unverified until it is. A parse failure on a real qptiff will correctly report `labels_unverifiable` rather than silently mislabeling — but "correctly refuses to guess" and "correctly reads real qptiff files" are different claims; only the first one has been checked.

## Format coverage (as of 2026-09-04)

| Format | Status |
|---|---|
| CODEX/PhenoCycler `.tif` (ImageJ `TCYX` hyperstack) | Fixed and verified — this is the real CRC TMA data. Pixel-level periodicity check passed by hand on one real core. |
| Plain flat `.tif` (`CYX`, already segmentation-ready) | Unaffected — `needs_flattening=False`, passes through unchanged. |
| `.ome.tif` / `.ome.tiff` | Name reading implemented and unit-tested on synthetic files (named, unnamed, 3-level pyramid). **Not yet run through the full segmentation pipeline on a real file** — `segmentation.py`'s glob still excludes `.ome.tif` by design (see `_find_masked_tifs`'s docstring); that exclusion has not been revisited as part of this change. |
| `.qptiff` | Name-reading code exists, **untested against a real file**. Pyramidal structure and RGB brightfield (`YXS` axes) are known additional wrinkles not yet handled — `inspect_stack()` will raise cleanly on `YXS` rather than mis-processing it, which is the correct behavior for now, not a fix. |

Extending real support to `.ome.tif`/`.qptiff` as first-class segmentation inputs (glob inclusion, pyramid-level selection, RGB handling) is intentionally deferred until after the first successful CODEX production run — see the 2026-09-04 conversation for the reasoning: validating the format this pipeline actually needs today took priority over generalizing to formats with no dataset behind them yet.

## Dependencies

`tifffile`, `numpy`. Both already required by `segmentation.py`; no new external dependency introduced.

## Notes / risks

- **Cycle-major flatten order is an assumption, not something read from metadata (confidence: high that it's correct for this dataset, verified by `verify_nuclear_periodicity`; would NOT be automatically caught if a future acquisition used the opposite convention).** `TCYX`/`ZCYX` flattening assumes the outer axis is cycle/z and reshapes as `(cycle, channel) → cycle*n_per_cycle + channel`. A hypothetical `CTYX` file (channel outer, cycle inner) has the identical shape tuple but the opposite memory layout — `inspect_stack()` does not special-case it and would raise `ValueError: Unhandled TIFF axis order`, which is the safe outcome, but if a future format silently reported itself as `TCYX` while actually being channel-major, only the periodicity check would catch it. That check is therefore load-bearing, not optional — don't run with `channel_check=False` on new acquisition types without independently confirming order first.
- **Panel identity is trusted, not verified, beyond periodicity (confidence: high, by design).** `verify_nuclear_periodicity` confirms the *nuclear* channels are self-consistent; it says nothing about whether channel 3 of cycle 7 is really `"CD45RA"` and not some other marker imaged in the wrong cycle. That level of verification requires either independent biological knowledge (does the CD45RA stain pattern look right?) or an acquisition log cross-check — outside what pixel statistics alone can confirm.
