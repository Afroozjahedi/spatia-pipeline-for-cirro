# QuPath: exporting ROI masks

**Where this fits:** a manual QuPath step, before ROI masking. Run this on annotations already drawn or detected in an open QuPath project — either by hand, or via the automatic TMA core-detection script below.

## What it does

Exports the mask images and label files that ROI masking needs: a mask image per slide, a text file listing every region with its position/size/area, and optionally individual cropped masks per region.

## How to run it

Open your QuPath project, paste the script into QuPath's Script Editor, and click Run. It prompts you for which annotations to include, resolution (leave at full resolution unless you have a specific reason not to), where to save the output, and whether to also save individual per-region masks.

## Things to know

- This has to be run by a person clicking through QuPath dialogs each time — there's no way to run it automatically or unattended.
- The area numbers this produces (used later for triad density) haven't yet been confirmed against a real QuPath run — worth a quick sanity check of the numbers and image names the first time you use this for a real density calculation.
- If you change the resolution setting away from full resolution, double-check with whoever runs ROI masking that coordinates still line up — at full resolution this isn't an issue.
