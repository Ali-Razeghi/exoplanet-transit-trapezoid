# transit_trapezoid_local.py
# --------------------------
# Educational, fast local transit light-curve generator (no internet required)
# Professional trapezoid model using orbital geometry.
#
# Inputs (expert-friendly):
#   - Orbital period P (days)
#   - Mid-transit time t0 (days; arbitrary zero-point)
#   - Rp/Rs (planet-to-star radius ratio)
#   - a/Rs (semi-major axis in stellar radii)
#   - Impact parameter b (0=centered, ~1=grazing)
#   - Total observing span (days), sampling interval (minutes), white noise (ppm)
#
# Model:
#   - Inclination from b = (a/Rs) * cos(i)  ->  sin(i) = sqrt(1 - (b/(a/Rs))^2)
#   - T14 = (P/pi) * asin[ (Rs/a) * sqrt((1+Rp/Rs)^2 - b^2) / sin(i) ]
#   - T23 = (P/pi) * asin[ (Rs/a) * sqrt((1-Rp/Rs)^2 - b^2) / sin(i) ]  (if argument valid; else no flat bottom)
#   - Trapezoid depth = (Rp/Rs)^2
#
# Plot:
#   - Single time-series light curve (normalized flux)
#   - Grid + tight layout
#
# Notes:
#   - If geometry is invalid (e.g., too large b), code warns and falls back to a simple box with 3% P.

import sys
import math
import numpy as np
import matplotlib
try:
    # Robust backend for PyCharm/Windows to avoid interagg tostring_rgb issues
    matplotlib.use("TkAgg")
except Exception:
    pass
import matplotlib.pyplot as plt

def ask_float(prompt, default=None, vmin=None, vmax=None):
    """Prompt a float with optional bounds and default."""
    txt = f"{prompt}"
    if default is not None:
        txt += f" [default {default}]"
    txt += ": "
    while True:
        s = input(txt).strip()
        if not s and default is not None:
            val = float(default)
        else:
            try:
                val = float(s)
            except Exception:
                print("Please enter a valid number.")
                continue
        if vmin is not None and val < vmin:
            print(f"Value must be >= {vmin}")
            continue
        if vmax is not None and val > vmax:
            print(f"Value must be <= {vmax}")
            continue
        return val

def trapezoid_model_time_series(t, P, t0, RpRs, aRs, b):
    """
    Build a trapezoid transit model using orbital geometry.
    Returns normalized flux array (1 outside transit).
    """
    # Validate geometry; compute sin(i)
    # b = (a/Rs) * cos(i)  ->  cos(i) = b / (a/Rs)
    small_eps = 1e-12
    if aRs <= 0:
        raise ValueError("a/Rs must be > 0.")
    cosi = b / aRs
    if abs(cosi) >= 1.0:
        # Impossible inclination; no transit. Return flat = 1.
        return np.ones_like(t, dtype=float), {"valid": False, "reason": "cos(i) >= 1 → no transit"}

    sini = math.sqrt(max(0.0, 1.0 - cosi * cosi))

    depth = RpRs * RpRs
    if depth <= 0:
        return np.ones_like(t, dtype=float), {"valid": False, "reason": "Rp/Rs <= 0 → no transit depth"}

    # Helper to compute duration term: (P/pi)*asin( (Rs/a) * sqrt( (1±RpRs)^2 - b^2 ) / sin(i) )
    def _dur_term(edge_factor):
        # edge_factor = 1 + RpRs (T14) or 1 - RpRs (T23)
        arg_sqrt = (edge_factor ** 2) - (b ** 2)
        if arg_sqrt <= 0:
            return None  # no solution (grazing beyond this limit)
        inside = (1.0 / aRs) * math.sqrt(arg_sqrt) / max(sini, small_eps)
        # Numerical safety: asin argument in [-1,1]
        inside = max(-1.0, min(1.0, inside))
        return (P / math.pi) * math.asin(inside)

    T14 = _dur_term(1.0 + RpRs)
    T23 = _dur_term(abs(1.0 - RpRs))  # abs for safety if RpRs > 1 (unphysical but robust)

    if T14 is None or T14 <= 0:
        # No geometric solution: return flat
        return np.ones_like(t, dtype=float), {"valid": False, "reason": "No geometric T14 (likely too large b)"}

    # If T23 invalid (e.g., grazing), set to 0 → pure triangle
    if (T23 is None) or (T23 < 0):
        T23 = 0.0

    # Ingress/egress duration:
    tau = 0.5 * (T14 - T23)
    tau = max(tau, 0.0)

    # Build trapezoid over ALL transits in the time window
    flux = np.ones_like(t, dtype=float)
    # transit centers within [t.min, t.max]
    n_min = math.floor((t[0] - t0) / P)
    n_max = math.ceil((t[-1] - t0) / P)
    centers = [t0 + n * P for n in range(n_min, n_max + 1)]

    for c in centers:
        # Regions relative to center
        x = t - c
        # |x| <= T23/2 → flat bottom (depth)
        if T23 > 0:
            in_flat = np.abs(x) <= (0.5 * T23)
            flux[in_flat] -= depth

        # Ingress: T23/2 < |x| <= (T23/2 + tau) → linearly down/up
        in_wings = (np.abs(x) > (0.5 * T23)) & (np.abs(x) <= (0.5 * T23 + tau))
        # Linear ramp height from 0 to depth across tau
        ramp = (np.abs(x[in_wings]) - 0.5 * T23) / max(tau, small_eps)  # 0..1
        flux[in_wings] -= depth * (1.0 - ramp)  # from depth to 0

        # Outside → no change (flux=1)

    meta = {
        "valid": True,
        "T14_days": T14,
        "T23_days": T23,
        "tau_days": tau,
        "depth": depth,
        "sini": sini,
    }
    return flux, meta

def main():
    print("=== Fast Transit Light Curve (Professional Trapezoid; No Internet) ===")

    # --- Expert inputs (with concise explanations) ---
    P   = ask_float("Orbital period P (days) (time to complete one orbit)", 3.0, vmin=0.01)
    t0  = ask_float("Mid-transit time t0 (days) (center time of a reference transit; arbitrary zero-point)", 1.5)
    RpRs = ask_float("Rp/Rs (planet radius / star radius)", 0.02, vmin=0.0001, vmax=0.5)
    aRs  = ask_float("a/Rs (semi-major axis / star radius)", 15.0, vmin=1.01)
    b    = ask_float("Impact parameter b (0=centered, ~1=grazing)", 0.2, vmin=0.0)

    total_days  = ask_float("Total time span (days) (length of simulated observation)", 20.0, vmin=0.1)
    cadence_min = ask_float("Sampling interval (minutes) (time step between data points)", 10.0, vmin=0.05)
    noise_ppm   = ask_float("Noise level (ppm) (random scatter added to the light curve)", 300.0, vmin=0.0)

    # --- Build time grid ---
    dt = float(cadence_min) / (60.0 * 24.0)  # days
    t = np.arange(0.0, float(total_days) + 1e-12, dt)

    # --- Generate trapezoid model ---
    flux, meta = trapezoid_model_time_series(t, P, t0, RpRs, aRs, b)

    # --- Add white noise if requested ---
    if noise_ppm > 0:
        flux = flux + np.random.normal(0.0, noise_ppm * 1e-6, size=flux.size)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(10, 4.6))
    ax.plot(t, flux, "-")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Normalized Flux")
    title = "Transit Light Curve (Trapezoid Model)"
    ax.set_title(title)
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

    # --- Print useful metrics for the user ---
    print("\n=== Derived transit geometry (trapezoid) ===")
    if meta.get("valid", False):
        print(f"T14 (first-to-fourth contact): {meta['T14_days']:.6f} days")
        print(f"T23 (second-to-third contact): {meta['T23_days']:.6f} days")
        print(f"Ingress/Egress duration tau:   {meta['tau_days']:.6f} days")
        print(f"Transit depth (Rp/Rs)^2:       {meta['depth']:.6f}")
        print(f"sin(i) from b, a/Rs:           {meta['sini']:.6f}")
    else:
        print(f"Model marked invalid: {meta.get('reason','unknown reason')}")
        print("Tip: ensure 0 ≤ b < a/Rs and Rp/Rs > 0.")

if __name__ == "__main__":
    main()
