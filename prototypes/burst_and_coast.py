import math
import os

import numpy as np

# Match launch/settingsBurstAndCoastFish.sh and BurstAndCoastFish.cpp
length = 0.2
amplitudeFactor = 1.0
T_bout = 0.5
lmbda = 0.2 * length
lmbda_bout = 2.0 * length
c = lmbda_bout / T_bout

a_x = [0, 0.2, 1]
a_y = [0.02, 0.01, 0.1]
c2, c1, c0 = np.polyfit(a_x, a_y, 2)


def a(x):
    n = x / length
    return amplitudeFactor * (c2 * n**2 + c1 * n + c0) * length


def da_dx(x):
    return amplitudeFactor * (2.0 * c2 * x / length + c1)


def _Y_pieces(x):
    x = np.asarray(x, dtype=float)
    scalar = x.ndim == 0
    xm = np.atleast_1d(x) % lmbda_bout
    xm = np.where(xm < 0, xm + lmbda_bout, xm)
    p_coast = lmbda_bout - lmbda
    p_rise = lmbda_bout - 0.75 * lmbda
    p_mid = lmbda_bout - 0.25 * lmbda
    m_rise = (xm >= p_coast) & (xm < p_rise)
    m_mid = (xm >= p_rise) & (xm < p_mid)
    m_fall = (xm >= p_mid) & (xm <= lmbda_bout)
    return scalar, xm, m_rise, m_mid, m_fall


def Y(x):
    scalar, xm, m_rise, m_mid, m_fall = _Y_pieces(x)
    out = np.zeros_like(xm)
    theta = (xm - lmbda_bout) / lmbda
    out[m_rise] = -0.5 * (1.0 - np.cos(4.0 * np.pi * theta[m_rise]))
    out[m_mid] = -np.sin(2.0 * np.pi * theta[m_mid])
    out[m_fall] = 0.5 * (1.0 - np.cos(4.0 * np.pi * theta[m_fall]))
    if scalar:
        return float(out[0])
    return out


def dY_dx(x):
    scalar, xm, m_rise, m_mid, m_fall = _Y_pieces(x)
    out = np.zeros_like(xm)
    theta = (xm - lmbda_bout) / lmbda
    k = 2.0 * np.pi / lmbda
    out[m_rise] = -k * np.sin(4.0 * np.pi * theta[m_rise])
    out[m_mid] = -k * np.cos(2.0 * np.pi * theta[m_mid])
    out[m_fall] = k * np.sin(4.0 * np.pi * theta[m_fall])
    if scalar:
        return float(out[0])
    return out


def y(x, t):
    return a(x) * Y(x - c * t)


def dy_dx(x, t):
    xi = x - c * t
    return da_dx(x) * Y(xi) + a(x) * dY_dx(xi)


def compute_integral_xend(t):
    cumulative = 0
    x_end = 0
    dx = 1e-4
    while cumulative < length:
        cumulative += math.sqrt(1 + (dy_dx(x_end, t)) ** 2) * dx
        x_end += dx
    return x_end


def make_video(out_path=None, n_periods=8, fps=30, x_max=None, slowdown=2.0):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation

    if out_path is None:
        out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "burst_and_coast.gif")
    if x_max is None:
        x_max = length

    n_x = 400
    x = np.linspace(0.0, x_max, n_x)
    t_end = n_periods * T_bout
    n_frames = int(round(fps * t_end * slowdown))
    t = np.linspace(0.0, t_end, n_frames, endpoint=False)
    env = a(x)
    y_lim = 1.15 * float(np.max(np.abs(env)))

    yt_grid = np.stack([y(x, ti) for ti in t], axis=0)

    fig, (ax_mid, ax_xt) = plt.subplots(
        2, 1, figsize=(8.0, 5.2), dpi=120, gridspec_kw={"height_ratios": [1.15, 1.0]}
    )
    fig.subplots_adjust(hspace=0.38, left=0.10, right=0.97, top=0.90, bottom=0.12)

    (line,) = ax_mid.plot(x, yt_grid[0], color="C0", lw=2.0)
    ax_mid.plot(x, env, color="0.55", ls="--", lw=0.9, label=r"$\pm a(x)$")
    ax_mid.plot(x, -env, color="0.55", ls="--", lw=0.9)
    ax_mid.axhline(0.0, color="0.8", lw=0.6)
    ax_mid.set_xlim(0.0, x_max)
    ax_mid.set_ylim(-y_lim, y_lim)
    ax_mid.set_xlabel(r"$x$")
    ax_mid.set_ylabel(r"$y(x,t)$")
    ax_mid.set_title(
        r"Burst-and-coast  $L={:g}$, $T_{{\mathrm{{bout}}}}={:g}$, $\lambda={:g}$, $\lambda_{{\mathrm{{bout}}}}={:g}$".format(
            length, T_bout, lmbda, lmbda_bout
        )
    )
    ax_mid.legend(loc="upper left", frameon=False, fontsize=9)

    im = ax_xt.imshow(
        yt_grid.T,
        origin="lower",
        aspect="auto",
        extent=[0.0, t_end, 0.0, x_max],
        cmap="RdBu_r",
        vmin=-y_lim,
        vmax=y_lim,
        interpolation="nearest",
    )
    t_cursor = ax_xt.axvline(0.0, color="k", lw=1.0)
    ax_xt.set_xlabel(r"$t$")
    ax_xt.set_ylabel(r"$x$")
    ax_xt.set_title(r"space–time  $y(x,t)$")
    fig.colorbar(im, ax=ax_xt, fraction=0.046, pad=0.02, label=r"$y$")
    time_text = ax_mid.text(
        0.98, 0.08, "", transform=ax_mid.transAxes, ha="right", va="bottom", fontsize=10
    )

    def update(k):
        line.set_ydata(yt_grid[k])
        t_cursor.set_xdata([t[k], t[k]])
        time_text.set_text(r"$t/T_{{\mathrm{{bout}}}}={:.2f}$".format(t[k] / T_bout))
        return line, t_cursor, time_text

    anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=int(1000 / fps), blit=True)
    anim.save(out_path, writer="pillow", fps=fps)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    path = make_video()
    print("saved", path)
