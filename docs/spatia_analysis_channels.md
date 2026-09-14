# Image format support (CODEX / OME-TIFF / QPTIFF)

Before segmentation runs, the pipeline needs to know how channels are laid out in your raw image files and match each one to a real marker name. This is handled automatically for the formats below — it's not a pipeline step you run yourself, just something segmentation relies on.

## Supported formats

- **CODEX/PhenoCycler raw hyperstacks (`.tif`)** — fully supported and validated against real CRC TMA data. The pipeline figures out the acquisition structure (cycles × channels) on its own and flattens it to one plane per marker.
- **Plain flat multi-channel `.tif`** (already one channel per marker) — passes through unchanged.
- **OME-TIFF (`.ome.tif` / `.ome.tiff`)** — channel names are read automatically from the file when present. Tested on sample files; not yet run through the full pipeline on a real OME-TIFF dataset, since none has come through yet.
- **QPTIFF (Akoya)** — channel-name reading is built, but hasn't been tried against a real `.qptiff` file. Treat it as unverified until it's actually run on one.

## What you still need to provide

A panel file (`channelnames.txt`) listing your marker names in stack order is required for every format — the pipeline checks this against the real file (channel count, plus the pixel-level check below) rather than guessing marker identity from the image itself, since that's genuinely not something that can be read off the pixels.

## The nuclear-stain consistency check

For CODEX-style data, the pipeline double-checks that the channel order it assumed is actually correct: the nuclear-stain channel from every imaging cycle should correlate with the others, since they're all re-imaging the same nuclei. If this check fails, the run stops rather than continuing on a possibly-scrambled channel order. This is the main safety net against a channel mixup that would otherwise be very hard to notice after the fact — don't turn it off (`channel_check: false`) on a new dataset unless you've independently confirmed the channel order some other way first.

## Things to know

- OME-TIFF and QPTIFF support exist so the pipeline isn't limited to one platform for the methods paper's sake — no dataset in either format has actually gone through segmentation yet. If a real file in either format comes in, budget some time to validate it the way CODEX was validated, before trusting the output.
- If a file's channel names can't be read from a format that normally carries them (a corrupted or unusual OME-TIFF/QPTIFF), the pipeline refuses to guess and stops with an error rather than silently applying your panel positionally. That's intentional — not something to work around.
