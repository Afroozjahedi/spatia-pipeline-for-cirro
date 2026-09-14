# Format conversion (raw scanner TIFFs → standard TIFF)

**Pipeline position:** Step 0 (first step, off by default). `raw QPTIFF/OME-TIFF` → **format conversion** → ROI masking → segmentation → ...

## What this step does

Converts vendor-specific TIFF formats (QPTIFF, OME-TIFF) straight off the slide scanner into a standard multi-page TIFF that the rest of the pipeline can read. The original file's metadata (channel names, etc.) is saved alongside each converted file as a JSON, since that information doesn't otherwise survive the format conversion. This step is only needed if your raw files are QPTIFF/OME-TIFF — see the image-format-support doc for what the pipeline already reads natively.

## Configuring it

- `paths.raw_image_dir` (required) — where your raw `.qptiff`/OME-TIFF files are
- `paths.converted_tif_dir` — where converted files go
- `tif_conversion.image_format` — file extension to convert (default `.qptiff`)
- `tif_conversion.preserve_metadata` — keep the metadata JSON (default on, recommended)
- `tif_conversion.overwrite` — redo files that were already converted (default off — files already done are skipped)

## What you'll get

- One `.tif` per raw file, plus a matching metadata JSON

## Things to know

- If your raw files are plain QPTIFF without embedded channel-name metadata, the metadata JSON may not include channel names for those files — the next step (ROI masking) already has a fallback for this, so it's handled, but don't be surprised if some of these JSONs look sparse.
- A file that fails to convert (corrupted, unreadable) is skipped and reported rather than stopping the whole batch — but the step only flags the case where *nothing at all* converted, not "most files failed." Worth a quick check that your output folder has roughly as many files as you expected.
