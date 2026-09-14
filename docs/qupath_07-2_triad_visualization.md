# QuPath: triad visualization

**Where this fits:** a manual, optional QuPath step for visually reviewing detected triads, after triad detection and the QuPath export step have both run.

**Status: not yet run against a real QuPath project.** Written and reviewed, but not tested on an actual slide yet — check the number of triads it renders against the export step's own printed count the first time you use it.

## What it does

Draws each member of a detected triad as a colored dot directly in your QuPath project — gold for the anchor cell, blue and red for the two partner cells — with the distances between them attached so you can inspect any triad by clicking on it. All three cells of one triad share a name so you can find them together.

## How to run it

Open the QuPath project for the specific image you want to review, edit the file path at the top of the script to point at that image's exported triad file, then run it from QuPath's Script Editor. Save the QuPath project afterward if you want to keep what was drawn.

## Things to know

- Make sure the file you point the script at actually matches the image you have open in QuPath — there's no automatic check for this, so pointing at the wrong file draws triads from a different image onto the wrong tissue with no warning.
- This hasn't been run on a real project yet, so treat the first use as a test and compare the count of rendered triads against what the export step reported.
