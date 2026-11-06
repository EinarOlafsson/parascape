# 92 · The Map module preview in 3D

**Status: open (requested 2026-11-06).** Einar: "the umap in map is 2d not 3d".

## The plan
The live preview projects to three dimensions and is drawn turning slowly, painted in software (no new GL
context: the window's GL map must not be disturbed, see the memory on sprite VBOs) with depth shading, so it
reads as the 3D map it previews. Drag rotates it.

## Gate
The preview is a rotatable 3D cloud that redraws as the settings change.
