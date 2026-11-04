# PROGRESS_playback (lane/playback)

Measured with web/scripts/playback-timing.mjs (Playwright, swiftshader, 21-frame default config, Coarse and Medium, Curvature and Length Fit).

Finding (before): the hypothesis "interpolation work exceeds the 60 ms tick" is NOT what happened. Display gaps were ~60-70 ms (tick-bound);
a synthesized frame cost ~5-12 ms. The unevenness came from the tick COUNT: n = ceil(jump/0.004) up to 40 ticks per true frame, so late
(large-jump) transitions took up to 720 ms while early ones took 60 ms. Smooth off was a flat 60 ms/frame, so smooth on was 3x slower overall.
Also found: finish() scheduled a 3.75 s setTimeout that called showEnd()/player.stop, cutting a user replay started right after a run (fixed).

Change: Player is rAF/wall-clock driven, FRAME_MS=450 per true frame for both smooth on/off; sub-frame count n = min(jump-based n, dur/(1.6*EMA render cost));
late renders skip sub-frames. Hermite plans cached per frame (WeakMap). Smooth off steps at the same boundaries.

| config | before mean/early/late/max (ms) | after mean/early/late/max (ms) |
|---|---|---|
| smooth on (all 4 cfgs) | 174 / 60 / 288 / 720 | 451-452 / 452-453 / 450 / 472-484 |
| smooth off | 60 / 60 / 60 / 64-80 | 451-452 / 452-453 / 450 / 466-476 |

Unresolved: sub-frames are few (38 over 20 transitions, same as before for most transitions; cap 12 for the big jumps vs 40 before), so large jumps look slightly less
smooth; headless swiftshader costs may differ from a real GPU. FRAME_MS is a constant (no speed control). No core/ changes.
