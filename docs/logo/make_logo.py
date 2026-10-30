"""GEM-PIC-RZ logo: writes the SVGs of docs/logo with the lettering converted to outlines.
Same layout, palette and fonts as the QUARZ logo (github.com/AlexMPukhov/QUARZ, docs/logo).
usage: FONTS=<dir with the two fonts> python3 make_logo.py <out dir>
PNGs: rendered from the SVGs with a headless browser at twice the size, e.g.
  msedge --headless --force-device-scale-factor=2 --window-size=820,320 --screenshot=out.png logo.svg"""
import math, os, sys
from fontTools.ttLib import TTFont
from fontTools.varLib.instancer import instantiateVariableFont
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
OUT = sys.argv[1]
F = os.environ.get('FONTS', 'fonts') + '/'   # SpaceGrotesk[wght].ttf, IBMPlexMono-Medium.ttf from github.com/google/fonts (OFL)
grot = instantiateVariableFont(TTFont(F + 'SpaceGrotesk[wght].ttf'), {"wght": 700})
mono = TTFont(F + 'IBMPlexMono-Medium.ttf')

def text_path(font, s, x, y, size, tracking_em):
    gs, cmap, upm = font.getGlyphSet(), font.getBestCmap(), font['head'].unitsPerEm
    sc = size / upm
    pen = SVGPathPen(gs)
    for ch in s:
        g = cmap[ord(ch)]
        gs[g].draw(TransformPen(pen, (sc, 0, 0, -sc, x, y)))
        x += gs[g].width * sc + tracking_em * size
    return pen.getCommands(), x

SP = "M0 -1 C 0.1 -0.1, 0.1 -0.1, 1 0 C 0.1 0.1, 0.1 0.1, 0 1 C -0.1 0.1, -0.1 0.1, -1 0 C -0.1 -0.1, -0.1 -0.1, 0 -1 Z"
ACC, AMBER = '#0E8FB8', '#F2A93B'

def poly(pts):
    return 'M' + ' L'.join(f'{x:.2f} {y:.2f}' for x, y in pts)

def pulse(xc, yc, w, amp, k, x0, x1, n=240):
    """laser pulse: carrier sin(k (x - xc)) under a Gaussian envelope; returns carrier and envelope paths"""
    xs = [x0 + (x1 - x0) * i / n for i in range(n + 1)]
    env = [(x, yc - amp * math.exp(-((x - xc) / w) ** 2)) for x in xs]
    car = [(x, yc + (y - yc) * math.sin(k * (x - xc))) for x, y in env]
    return poly(car), poly(env + [(x, 2 * yc - y) for x, y in reversed(env)]) + ' Z'

def mark(ink):
    """letter G as a ring (gap from 8 to 52 deg) with its crossbar on the axis; inside: laser pulse,
    nested bubble sheaths of the wake behind it, witness bunch in the first bubble"""
    s = f'<g fill="none" stroke="{ink}">'
    for a, w, o in [(16, 1.8, .75), (23, 2.3, .62), (30, 2.8, .52), (37, 3.3, .44)]:
        cx = 104 - a    # every sheath closes at the pulse, x = 104
        s += (f'<path d="M104 100 C {104 - .3 * a:.1f} {100 - a:.1f}, {cx - .55 * a:.1f} {100 - a:.1f}, {cx - a:.1f} 100 '
              f'C {cx - .55 * a:.1f} {100 + a:.1f}, {104 - .3 * a:.1f} {100 + a:.1f}, 104 100 Z" stroke-width="{w}" opacity="{o}"/>')
    s += '</g>'
    r, (g0, g1) = 84, (math.radians(52), math.radians(368))
    s += (f'<path d="M{100 + r * math.cos(g0):.2f} {100 - r * math.sin(g0):.2f} A{r} {r} 0 1 0 {100 + r * math.cos(g1):.2f} '
          f'{100 - r * math.sin(g1):.2f}" fill="none" stroke="{ink}" stroke-width="14"/>')
    s += f'<rect x="150" y="93" width="41" height="14" fill="{ink}"/>'
    car, env = pulse(127, 100, 11, 20, 0.75, 104, 150)
    s += f'<path d="{env}" fill="{ACC}" opacity="0.22"/><path d="{car}" fill="none" stroke="{ACC}" stroke-width="4.2" stroke-linecap="round"/>'
    s += ''.join(f'<ellipse cx="{c}" cy="100" rx="{a}" ry="{b}" fill="{ACC}"/>' for c, a, b in [(76, 4.6, 6.5), (87, 4.2, 5.6)])
    return s

def svg(w, h, body, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-label="{title}">'
            f'<title>{title}</title>{body}</svg>\n')

def logo(dark):
    bg, ink, sub = ('#0B0F18', '#EEF1F4', '#A9B2BE') if dark else ('#F4F5F2', '#1F3A5F', '#566E8E')
    W, H = 820, 320
    b = f'<rect width="{W}" height="{H}" fill="{bg}"/>'
    b += f'<g><path d="{SP}" transform="translate(690 100) scale(150)" fill="{AMBER}" opacity="0.22"/>'
    b += f'<path d="{SP}" transform="translate(768 280) scale(22)" fill="{AMBER}" opacity="0.85"/>'
    b += f'<path d="{SP}" transform="translate(560 26) scale(11)" fill="{AMBER}" opacity="0.7"/></g>'
    b += f'<svg x="64" y="60" width="200" height="200" viewBox="-6 -6 212 212">{mark(ink)}</svg>'
    word, _ = text_path(grot, 'GEM-PIC-RZ', 300, 174, 76, 0.02)
    tag, _ = text_path(mono, 'GALILEAN ELECTROMAGNETIC PIC IN RZ', 302, 210, 13, 0.14)
    b += f'<path d="{word}" fill="{ink}"/><path d="{tag}" fill="{sub}"/>'
    return svg(W, H, b, 'GEM-PIC-RZ – Galilean electromagnetic PIC in RZ')

files = {'gempic-rz-logo-dark.svg': logo(True), 'gempic-rz-logo-light.svg': logo(False)}
for n, s in files.items():
    open(os.path.join(OUT, n), 'w', encoding='utf-8').write(s)
print('\n'.join(files))
