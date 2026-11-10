"""Exact Ellis null-ray map, illustrative skies. Author: Biswajit Jana.

This generates a scientific visualization from equations, not an image edit.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import ellipk, ellipkinc
from scipy.integrate import quad
from scipy.ndimage import gaussian_filter, map_coordinates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT/"figures/study-003"
OUT = ROOT/"results/study-003"


def ray_angle(q, observer=8.):
    """Total polar advance for inward rays; elliptic functions use parameter m."""
    q = np.asarray(q, dtype=float)
    if np.any(q < 0) or np.any(q >= np.sqrt(1+observer**2)) or np.any(q == 1):
        raise ValueError("Require inward, noncritical rays within observer's local light cone")
    result = np.empty_like(q)
    through = q < 1
    v = q[through]
    result[through] = v*(2*ellipk(v*v)-ellipkinc(np.arctan(1/observer), v*v))
    v = q[~through]
    result[~through] = 2*ellipk(1/(v*v))-ellipkinc(np.arcsin(v/np.sqrt(1+observer**2)), 1/(v*v))
    return result


def sky(seed, opposite=False):
    """Procedural RGB radiance: finite-size stars and diffuse bands, no catalog."""
    rng = np.random.default_rng(seed)
    h, w = 1024, 2048
    lat = np.linspace(np.pi/2, -np.pi/2, h)[:, None]
    lon = np.linspace(-np.pi, np.pi, w, endpoint=False)[None, :]
    bandlat = lat-0.45*np.sin(lon+(1.2 if opposite else 0.2))
    band = np.exp(-(bandlat/.12)**2)
    clouds = np.zeros((h, w))
    for scale, weight in [(3, .2), (12, .3), (40, .5)]:
        field = gaussian_filter(rng.normal(size=(h, w)), scale, mode="wrap")
        clouds += weight*field/field.std()
    dust = np.exp(-.8*np.maximum(clouds, 0))
    diffuse = .045*band*dust*(1+.25*np.tanh(clouds))
    color = np.array([.65, .8, 1.]) if opposite else np.array([1., .8, .58])
    rgb = .0005+diffuse[:, :, None]*color
    stars = np.zeros((h, w, 3))
    count = 35000
    latitude = np.arcsin(rng.uniform(-1, 1, count))
    longitude = rng.uniform(-np.pi, np.pi, count)
    rows = np.minimum(h-1, ((np.pi/2-latitude)/np.pi*(h-1)).astype(int))
    cols = ((longitude+np.pi)/(2*np.pi)*w).astype(int)%w
    brightness = np.minimum(rng.pareto(1.6, count)+.1, 30)*.06
    tint = rng.uniform(0, 1, count)
    colors = np.column_stack([.65+.35*tint, .75+.15*tint, 1-.4*tint])
    for channel in range(3):
        np.add.at(stars[:, :, channel], (rows, cols), brightness*colors[:, channel])
        stars[:, :, channel] = gaussian_filter(stars[:, :, channel], .65, mode="wrap")
    return rgb+stars


def sample_sky(texture, advance, azimuth):
    x = np.sin(advance)*np.cos(azimuth)
    y = np.sin(advance)*np.sin(azimuth)
    z = np.cos(advance)
    lon, lat = np.arctan2(y, x), np.arcsin(np.clip(z, -1, 1))
    h, w, _ = texture.shape
    coords = [(np.pi/2-lat)/np.pi*(h-1), (lon+np.pi)/(2*np.pi)*w]
    # Latitude is clamped by the coordinates; longitude wraps at its seam.
    return np.stack([map_coordinates(texture[:, :, c], coords, order=1, mode="grid-wrap") for c in range(3)], axis=-1)


def camera(width, height, samples, sky_a, sky_b, fov=38., observer=8.):
    radius = np.sqrt(1+observer**2)
    halfwidth = np.tan(np.deg2rad(fov/2))
    halfheight = halfwidth*height/width
    result = np.zeros((height, width, 3))
    for sy in range(samples):
        for sx in range(samples):
            x = ((np.arange(width)+(sx+.5)/samples)/width*2-1)*halfwidth
            y = (1-(np.arange(height)+(sy+.5)/samples)/height*2)*halfheight
            xx, yy = np.meshgrid(x, y)
            psi = np.arctan(np.hypot(xx, yy))
            q = radius*np.sin(psi)
            az = np.arctan2(yy, xx)
            advance = ray_angle(q, observer)
            through = q < 1
            for mask, texture in [(~through, sky_a), (through, sky_b)]:
                result[mask] += sample_sky(texture, advance[mask], az[mask])
    return result/samples**2


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    validations = []
    obs = 8.
    for q in [0., .1, .5, .9, .99, 1.01, 1.1, 2., 4.]:
        if q < 1:
            fn = lambda z: q/np.sqrt(1-q*q*np.cos(z)**2)
            direct = quad(fn, -np.pi/2, np.arctan(obs), epsabs=1e-10, epsrel=1e-10)[0]
        else:
            # Independent substitution x=sqrt(q²-1)/cos(t), regular at turning point.
            turn = np.sqrt(q*q-1)
            fn = lambda t: q/np.sqrt(turn*turn+np.cos(t)**2)
            direct = quad(fn, 0, np.pi/2, epsabs=1e-10)[0]+quad(fn, 0, np.arccos(turn/obs), epsabs=1e-10)[0]
        exact = float(ray_angle(q, obs))
        validations.append(dict(b_over_a=q, elliptic_angle=exact, quadrature_angle=direct, absolute_error=abs(exact-direct)))
    assert max(r["absolute_error"] for r in validations) < 1e-8
    # For a constant radiance on both skies, every ray returns that radiance.
    flat = np.ones((20, 40, 3))*.4
    uniform = camera(40, 28, 1, flat, flat)
    assert np.max(abs(uniform-.4)) < 1e-12
    a, b = sky(20261113), sky(20261114, opposite=True)
    radiance = camera(1200, 800, 2, a, b)
    display = np.clip(1-np.exp(-radiance*12), 0, 1)**(1/2.2)
    plt.imsave(FIG/"01-ideal-camera.png", display)
    fig, ax = plt.subplots(figsize=(12, 8), facecolor="#040811")
    ax.imshow(display)
    ax.axis("off")
    ax.set_title("ELLIS WORMHOLE · hypothetical camera view", color="white", fontsize=18, pad=14)
    label = dict(color="white", fontsize=11, bbox=dict(facecolor="#08101b", alpha=.85, edgecolor="none", pad=5))
    arrow = dict(arrowstyle="-", color="#f4d59f", lw=1.4)
    critical = np.rad2deg(np.arcsin(1/np.sqrt(65)))
    pixel_radius = 600*np.tan(np.deg2rad(critical))/np.tan(np.deg2rad(19))
    ax.annotate("A · Light from the other end\nColor and star field are illustrative", xy=(605, 420), xytext=(760, 645), arrowprops=arrow, **label)
    ax.annotate("B · Critical curve / winding images\nb = a; angular radius ≈ 7.13°", xy=(600+pixel_radius, 400), xytext=(790, 105), arrowprops=arrow, **label)
    ax.annotate("C · Our-side sky, gravitationally remapped", xy=(220, 300), xytext=(35, 70), arrowprops=arrow, **label)
    fig.text(.5, .04, "Exact null-ray mapping · observer x/a = 8 · horizontal field 38°\nSynthetic skies · ideal camera · no accretion, detector model or observational claim", ha="center", color="#bfcad9", fontsize=11)
    fig.tight_layout(rect=[0, .075, 1, 1])
    fig.savefig(FIG/"02-camera-with-legends.png", dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)
    with (OUT/"geodesic-validation.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(validations[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(validations)
    qgrid = np.r_[np.linspace(0, .99, 100), np.linspace(1.01, 4., 100)]
    with (OUT/"ray-map.csv").open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["b_over_a", "advance_radians", "destination"])
        writer.writerows((q, float(ray_angle(q)), "other_end" if q < 1 else "observer_side") for q in qgrid)
    config = {"author": "Biswajit Jana", "model": "ultrastatic symmetric Ellis null geodesics",
        "observer_x_over_a": 8, "horizontal_fov_degrees": 38, "pixels": [1200, 800], "subpixels_per_pixel": 4,
        "sky_seeds": [20261113, 20261114], "sky_source": "procedural finite-size stars and diffuse bands; no catalog",
        "critical_angular_radius_degrees": float(critical), "display": "[1-exp(-12 radiance)]^(1/2.2)",
        "uniform_radiance_max_error": float(np.max(abs(uniform-.4))),
        "max_angle_quadrature_error": max(r["absolute_error"] for r in validations),
        "limitations": ["finite sampling does not resolve arbitrarily high winding images", "no detector noise or PSF",
                        "no radiation transport through matter", "no realistic source or observer distance inferred",
                        "surface brightness conserved; no artificial magnification multiplier"]}
    (OUT/"summary.json").write_text(json.dumps(config, indent=2)+"\n", newline="\n")
    paths = [p for folder in (OUT, FIG) for p in sorted(folder.iterdir()) if p.suffix in (".json", ".csv", ".png") and p.name != "sha256.json"]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (OUT/"sha256.json").write_text(json.dumps(hashes, indent=2)+"\n", newline="\n")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
