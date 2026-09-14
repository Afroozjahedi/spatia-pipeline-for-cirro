# ROI masking

**Pipeline position:** Step 1 (off by default). Format conversion → **ROI masking** → segmentation → preprocessing → ...

## What this step does

Crops out the regions of interest (ROIs) you defined in QuPath — either manually drawn ROIs or an automatic TMA core detection — from your whole-slide/whole-core images, using the masks QuPath exported. Produces one cropped, masked image per ROI, ready for segmentation.

## Configuring it

- `experiment.groups` — reused from elsewhere in your config
- `paths.converted_tif_dir` — input images
- `paths.qupath_roi_dir` (required) — the folder QuPath exported your masks and ROI labels into
- `paths.masked_roi_dir` — where the cropped ROI images go

## What you'll get

- One masked, cropped image per ROI (plus a preview PNG and a channel-info text file for each)
- `roi_processing_summary.csv` summarizing what was processed

## Things to know

- If a single channel fails to read while cropping a ROI, that channel is left blank (all zeros) in the output rather than failing the whole ROI — a warning prints when this happens, so watch the run's console output rather than assuming a blank channel later on means "no real signal."
- Very large whole-slide images are allowed through without a file-size safety check (this is necessary, since real whole-slide images are legitimately huge) — a corrupted or unexpectedly oversized file may use a lot of memory before anything else catches the problem, rather than failing right away.
- If your filenames don't clearly match your `experiment.groups` prefixes, ROI matching can behave unexpectedly — worth a quick visual check of the first run's output against your filename conventions if any slide IDs look off.
