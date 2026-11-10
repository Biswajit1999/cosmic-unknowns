# Study 003: hypothetical camera data

**Author: Biswajit Jana**

This is an exact null-direction map in the prescribed symmetric ultrastatic Ellis geometry, rendered against procedural RGB source skies. There is no observed wormhole image or astronomical source catalog in this dataset.

- `ray-map.csv`: 200 impact parameters b/a, total angular advances in radians, and source-side classifications; the critical value b/a=1 is excluded.
- `geodesic-validation.csv`: nine direct-quadrature comparisons across crossing and returning rays. Values are deterministic; statistical error bars would be inappropriate.
- `summary.json`: camera and source parameters, seeds, tone mapping, angular radius, checks and limitations.
- `sha256.json`: hashes of data and the two calculated image files.

Run `scripts/render_wormhole_camera.py`. The unlabeled image is `figures/study-003/01-ideal-camera.png`; its explanatory companion is `figures/study-003/02-camera-with-legends.png`. The [report](../../reports/006-what-a-wormhole-could-look-like.md) derives the map and distinguishes predicted ray directions from assumed source colors and brightnesses.

No image-generation model modifies either calculated image. AI-generated photorealistic interpretations, if present, are stored separately under `figures/concept/` and do not constitute quantitative data.
