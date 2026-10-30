# GEM-PIC-RZ logo

A laser wakefield on the axis of the code's quasi-cylindrical (r, z) geometry. The ring is the G; its
crossbar is the axis. On the axis, the laser pulse (carrier under a Gaussian envelope) drives the nested
bubble sheaths of its wake, and a witness bunch sits in the first bubble. The style follows the
[QUARZ](https://github.com/AlexMPukhov/QUARZ) logo. The amber sparks in the background mark that the
code was developed with an AI assistant (Claude).

| file | use |
|---|---|
| `gempic-rz-logo-dark.svg` / `.png` | full logo on a dark background (README in dark mode, slides) |
| `gempic-rz-logo-light.svg` / `.png` | full logo on a light background (README, posters, papers) |
| `gempic-rz-mark-dark.svg` / `-light.svg` (`.png`) | the G alone, transparent background |
| `gempic-rz-avatar-dark.svg` / `.png` | 512 × 512, bolder mark: GitHub avatar, social media |
| `gempic-rz-avatar-light.svg` / `.png` | the same on the light background, navy ring |

Colours: background `#0B0F18` (dark) / `#F4F5F2` (light), ink `#EEF1F4` / navy `#1F3A5F` (tagline `#566E8E`
light, `#A9B2BE` dark), laser and bunch cyan `#0E8FB8`, sparks amber `#F2A93B`. Lettering: Space Grotesk Bold
(wordmark) and IBM Plex Mono Medium (tagline), both SIL Open Font License, converted to outlines, so the SVGs
need no fonts. `make_logo.py` regenerates the SVGs (`FONTS=<dir> python3 make_logo.py <out dir>`); the PNGs are
rendered from them with a headless browser (logo and mark at 2×, the marks on a transparent background).
