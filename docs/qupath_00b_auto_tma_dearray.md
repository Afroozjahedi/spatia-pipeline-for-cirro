# QuPath: automatic TMA core detection

**Where this fits:** a manual QuPath step, before ROI masking — an automatic alternative to drawing TMA core annotations by hand. Run this first, then run the ROI-mask-export script above on its output. TMA cohorts only, not for whole-slide images.

**Status: not yet run against a real slide.** This script is written and ready to try, but hasn't been validated on an actual scan yet — treat its first real run as a test, and visually check the detected core grid against the slide before trusting it.

## What it does

Automatically detects TMA cores on a whole-array scan using QuPath's built-in core-detection tool, instead of drawing each core boundary by hand.

## How to run it

From the command line (QuPath supports running scripts without opening the GUI) or from inside QuPath's Script Editor. Before running, edit a few settings directly in the script:

- the core diameter, from your TMA's construction specs — don't guess this, a wrong value will misdetect cores
- the grid size and labeling order (how many rows/columns, and whether labels run row-first or column-first)
- a density threshold that affects how strictly it distinguishes real tissue from empty background — the default is the same value QuPath's own examples use, but it hasn't been tuned or checked against these slides yet

## Things to know

- Find the right input file first. Some of the per-core files already saved from earlier work on this project look like they were exported *after* automatic dearraying already happened, not the original whole-array scan this script needs — worth locating (or re-exporting) the actual pre-dearray scan before the first real test.
- Since this hasn't been run yet, plan to visually spot-check the detected grid against the real slide the first time you use it, rather than trusting the output blindly.
