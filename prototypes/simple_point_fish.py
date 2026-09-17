"""
Wang et al. 2022 burst-and-coast fish school
(PLoS Comput Biol 18:e1009437), plus a StefanFish-like PID prototype.

Wang: at each kick, δφ = Σ_k MINs (Att+Ali) + γ_R g, then coast
       v(t) = v0 exp(-t/τ0).

StefanPID: same Wang social setpoints on a kick clock, but continuous
           planar act()-style controls tracked by PD:
             b = turn bias (dimless, >0 → CCW)
             a = period factor (dimless, T=T0(1+a); a<0 → faster)

Fig. 5A schooling: (γ_Att, γ_Ali) = (0.04, 0.20)
Fig. 5B milling:   (γ_Att, γ_Ali) = (0.05, 0.05)
"""

from __future__ import print_function

import math
import os
import time as _time

import numpy as np


def _wrap(a):
    return ((a + math.pi) % (2.0 * math.pi)) - math.pi


def _wrap_np(a):
    return ((a + np.pi) % (2.0 * np.pi)) - np.pi


def _wang_delta_phi_all(
    i, x, y, phi, gamma_att, gamma_ali, d_att, d_ali, l_att, l_ali
):
    """Pairwise δφ from fish i (Wang eqs 7–8, 11–16)."""
    dx = x - x[i]
    dy = y - y[i]
    d = np.maximum(np.hypot(dx, dy), 1e-6)
    psi = _wrap_np(np.arctan2(dy, dx) - phi[i])
    dphi = _wrap_np(phi - phi[i])
    f_att = gamma_att * (d / d_att - 1.0) / (1.0 + (d / l_att) ** 2)
    f_ali = gamma_ali * (d / d_ali + 1.0) * np.exp(-((d / l_ali) ** 2))
    o_att = 1.395 * np.sin(psi) * (1.0 - 0.33 * np.cos(psi))
    e_att = 0.9326 * (1.0 - 0.48 * np.cos(dphi) - 0.31 * np.cos(2.0 * dphi))
    e_ali = 0.9012 * (1.0 + 0.6 * np.cos(psi) - 0.32 * np.cos(2.0 * psi))
    o_ali = 1.6385 * np.sin(dphi) * (1.0 + 0.3 * np.cos(2.0 * dphi))
    out = f_att * o_att * e_att + f_ali * e_ali * o_ali
    out[i] = 0.0
    return out


def _disk_ic(n, l_att, rng):
    R = 0.5 * l_att * math.sqrt(n / math.pi)
    rad = np.sqrt(rng.random(n)) * R
    th = 2.0 * math.pi * rng.random(n)
    return rad * np.cos(th), rad * np.sin(th), 2.0 * math.pi * rng.random(n) - math.pi


class WangBurstCoastSimulator(object):
    """Open-world burst-and-coast school (Wang / Calovi / Lei)."""

    def __init__(
        self,
        n=100,
        k=1,
        gamma_att=0.04,
        gamma_ali=0.20,
        l_att_cm=28.0,
        l_ali_cm=28.0,
        d_att_cm=3.0,
        d_ali_cm=6.0,
        gamma_r=0.20,
        v0_cm=14.0,
        tau0_s=0.8,
        kick_len_cm=7.0,
        seed=71,
    ):
        self.n = int(n)
        self.k = max(1, int(k))
        self.gamma_att = float(gamma_att)
        self.gamma_ali = float(gamma_ali)
        self.l_att = float(l_att_cm)
        self.l_ali = float(l_ali_cm)
        self.d_att = float(d_att_cm)
        self.d_ali = float(d_ali_cm)
        self.gamma_r = float(gamma_r)
        self.v0_cm = float(v0_cm)
        self.tau0 = float(tau0_s)
        self.kick_len_cm = float(kick_len_cm)
        self.time = 0.0
        self.rng = np.random.RandomState(int(seed))

        self.x, self.y, self.phi = _disk_ic(self.n, self.l_att, self.rng)
        self.speed = np.full(self.n, self.v0_cm)
        self.t_kick = np.zeros(self.n)
        self.tau_kick = np.full(self.n, 0.5)
        self.v_kick = np.full(self.n, self.v0_cm)
        for i in range(self.n):
            self._sample_kick(i, defer=True)
            self.t_kick[i] = -self.rng.uniform(0.0, 0.9 * self.tau_kick[i])
            age = self.time - self.t_kick[i]
            self.speed[i] = self.v_kick[i] * math.exp(-age / self.tau0)

    def _delta_phi_all(self, i):
        return _wang_delta_phi_all(
            i, self.x, self.y, self.phi,
            self.gamma_att, self.gamma_ali, self.d_att, self.d_ali, self.l_att, self.l_ali,
        )

    def _sample_kick(self, i, defer=False):
        v0 = float(np.clip(self.rng.normal(self.v0_cm, 0.25 * self.v0_cm), 4.0, 28.0))
        length = float(max(2.0, self.rng.normal(self.kick_len_cm, 0.25 * self.kick_len_cm)))
        length = min(length, 0.98 * v0 * self.tau0)
        tau = -self.tau0 * math.log(max(1e-6, 1.0 - length / (v0 * self.tau0)))
        self.tau_kick[i] = max(0.15, min(tau, 2.5))
        self.v_kick[i] = v0
        if not defer:
            self.speed[i] = v0
            self.t_kick[i] = self.time

    def _do_kick(self, i):
        dphi = self._delta_phi_all(i)
        if self.k >= self.n - 1:
            social = float(dphi.sum())
        else:
            idx = np.argpartition(np.abs(dphi), -self.k)[-self.k :]
            social = float(dphi[idx].sum())
        self.phi[i] = _wrap(float(self.phi[i]) + social + self.gamma_r * float(self.rng.normal()))
        self._sample_kick(i)

    def step(self, dt):
        dt = float(dt)
        age = self.time - self.t_kick
        self.speed = self.v_kick * np.exp(-age / self.tau0)
        self.x = self.x + self.speed * np.cos(self.phi) * dt
        self.y = self.y + self.speed * np.sin(self.phi) * dt
        for i in np.nonzero(age + dt >= self.tau_kick)[0]:
            self._do_kick(int(i))
        self.time += dt

    def run(self, dt, n_steps):
        for _ in range(int(n_steps)):
            self.step(dt)

    def polarization(self):
        return float(np.hypot(np.mean(np.cos(self.phi)), np.mean(np.sin(self.phi))))

    def milling(self):
        vx = self.speed * np.cos(self.phi)
        vy = self.speed * np.sin(self.phi)
        theta = np.arctan2(self.y - self.y.mean(), self.x - self.x.mean())
        phi_bar = np.arctan2(vy - vy.mean(), vx - vx.mean())
        return float(np.abs(np.mean(np.sin(_wrap_np(phi_bar - theta)))))


class StefanPIDSchool(object):
    """
    Prototype of CUP3D planar act() + outer PD.

    act ~ {b, a}:
      b = turn/curvature bias (dimless; >0 CCW)
      a = period factor       (dimless; T=T0(1+a); a<0 → faster)

    Wang social refreshes (φ_des, v_des) ~every 0.5 s; PD tracks them smoothly.
    """

    def __init__(
        self,
        n=25,
        k=1,
        gamma_att=0.04,
        gamma_ali=0.20,
        l_att_cm=28.0,
        l_ali_cm=28.0,
        d_att_cm=3.0,
        d_ali_cm=6.0,
        gamma_r=0.20,
        v0_cm=14.0,
        decide_s=0.5,
        yaw_rate_gain=2.5,
        speed_tau_s=0.35,
        kp_phi=1.2,
        kd_phi=0.25,
        kp_v=0.08,
        b_max=1.5,
        a_max=0.5,
        seed=71,
    ):
        self.n = int(n)
        self.k = max(1, int(k))
        self.gamma_att = float(gamma_att)
        self.gamma_ali = float(gamma_ali)
        self.l_att = float(l_att_cm)
        self.l_ali = float(l_ali_cm)
        self.d_att = float(d_att_cm)
        self.d_ali = float(d_ali_cm)
        self.gamma_r = float(gamma_r)
        self.v0_cm = float(v0_cm)
        self.decide_s = float(decide_s)
        self.yaw_rate_gain = float(yaw_rate_gain)
        self.speed_tau_s = float(speed_tau_s)
        self.kp_phi = float(kp_phi)
        self.kd_phi = float(kd_phi)
        self.kp_v = float(kp_v)
        self.b_max = float(b_max)
        self.a_max = float(a_max)
        self.time = 0.0
        self.rng = np.random.RandomState(int(seed))

        self.x, self.y, self.phi = _disk_ic(self.n, self.l_att, self.rng)
        self.speed = np.full(self.n, self.v0_cm)
        self.omega = np.zeros(self.n)
        self.phi_des = self.phi.copy()
        self.v_des = np.full(self.n, self.v0_cm)
        self.t_decide = -self.rng.uniform(0.0, self.decide_s, size=self.n)
        self.tau_decide = np.full(self.n, self.decide_s)
        self.b = np.zeros(self.n)
        self.a = np.zeros(self.n)

    def _social(self, i):
        dphi = _wang_delta_phi_all(
            i, self.x, self.y, self.phi,
            self.gamma_att, self.gamma_ali, self.d_att, self.d_ali, self.l_att, self.l_ali,
        )
        if self.k >= self.n - 1:
            social = float(dphi.sum())
        else:
            idx = np.argpartition(np.abs(dphi), -self.k)[-self.k :]
            social = float(dphi[idx].sum())
        return social + self.gamma_r * float(self.rng.normal())

    def _refresh_setpoint(self, i):
        self.phi_des[i] = _wrap(float(self.phi[i]) + self._social(i))
        self.v_des[i] = float(np.clip(self.rng.normal(self.v0_cm, 0.25 * self.v0_cm), 6.0, 22.0))
        self.tau_decide[i] = float(np.clip(self.rng.normal(self.decide_s, 0.12), 0.25, 1.2))
        self.t_decide[i] = self.time

    def _pd_act(self, i):
        e_phi = _wrap(float(self.phi_des[i] - self.phi[i]))
        b = self.kp_phi * e_phi - self.kd_phi * float(self.omega[i])
        self.b[i] = float(np.clip(b, -self.b_max, self.b_max))
        e_v = float(self.v_des[i] - self.speed[i])
        # a<0 → shorter period → higher speed in the plant
        self.a[i] = float(np.clip(-self.kp_v * e_v, -self.a_max, self.a_max))

    def step(self, dt):
        dt = float(dt)
        for i in range(self.n):
            if self.time - self.t_decide[i] >= self.tau_decide[i]:
                self._refresh_setpoint(i)
            self._pd_act(i)
            # Plant surrogate of act({b, a})
            self.omega[i] = self.yaw_rate_gain * self.b[i]
            self.phi[i] = _wrap(float(self.phi[i]) + self.omega[i] * dt)
            v_tgt = self.v0_cm / max(1.0 + self.a[i], 0.35)
            self.speed[i] += (v_tgt - self.speed[i]) * (dt / self.speed_tau_s)
            self.x[i] += self.speed[i] * math.cos(self.phi[i]) * dt
            self.y[i] += self.speed[i] * math.sin(self.phi[i]) * dt
        self.time += dt

    def run(self, dt, n_steps):
        for _ in range(int(n_steps)):
            self.step(dt)

    def polarization(self):
        return float(np.hypot(np.mean(np.cos(self.phi)), np.mean(np.sin(self.phi))))

    def milling(self):
        vx = self.speed * np.cos(self.phi)
        vy = self.speed * np.sin(self.phi)
        theta = np.arctan2(self.y - self.y.mean(), self.x - self.x.mean())
        phi_bar = np.arctan2(vy - vy.mean(), vx - vx.mean())
        return float(np.abs(np.mean(np.sin(_wrap_np(phi_bar - theta)))))


def _render_gif(out_path, xs, ys, phis, dt, burn_in_s, fps, title_prefix):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation, cm
    from matplotlib.patches import Polygon

    n_rec, n_fish = xs.shape
    P_t = np.hypot(np.mean(np.cos(phis), axis=1), np.mean(np.sin(phis), axis=1))
    theta = np.arctan2(ys - ys.mean(1, keepdims=True), xs - xs.mean(1, keepdims=True))
    M_t = np.abs(np.mean(np.sin(_wrap_np(phis - theta)), axis=1))
    print("  P={:.3f}  M={:.3f}".format(P_t.mean(), M_t.mean()))

    stride = max(1, int(round(1.0 / (fps * dt))))
    frame_ids = list(range(0, n_rec, stride))
    trail = max(1, int(round(1.0 / dt)))
    colors = [cm.get_cmap("hsv", n_fish)(i) for i in range(n_fish)]

    fig, ax = plt.subplots(figsize=(7, 7), dpi=100)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    nose, half_w = 2.2, 0.7
    trails = [ax.plot([], [], "-", lw=0.7, alpha=0.4, color=colors[i])[0] for i in range(n_fish)]
    bodies = []
    for i in range(n_fish):
        p = Polygon(
            [[0, 0], [0, 0], [0, 0]], closed=True,
            facecolor=colors[i], edgecolor="0.15", lw=0.25, zorder=5,
        )
        ax.add_patch(p)
        bodies.append(p)
    title = ax.set_title("")

    def wedge(x, y, phi):
        c, s = math.cos(phi), math.sin(phi)
        tip = (x + nose * c, y + nose * s)
        left = (x - 0.4 * nose * c - half_w * s, y - 0.4 * nose * s + half_w * c)
        right = (x - 0.4 * nose * c + half_w * s, y - 0.4 * nose * s - half_w * c)
        return [tip, left, right]

    def view(fi):
        cx, cy = float(xs[fi].mean()), float(ys[fi].mean())
        half = 0.65 * max(float(np.ptp(xs[fi])), float(np.ptp(ys[fi])), 20.0) + 8.0
        ax.set_xlim(cx - half, cx + half)
        ax.set_ylim(cy - half, cy + half)

    def update(k):
        fi = frame_ids[k]
        i0 = max(0, fi + 1 - trail)
        for i in range(n_fish):
            trails[i].set_data(xs[i0 : fi + 1, i], ys[i0 : fi + 1, i])
            bodies[i].set_xy(wedge(float(xs[fi, i]), float(ys[fi, i]), float(phis[fi, i])))
        view(fi)
        title.set_text(
            "{}  N={}  t={:.0f}s  P={:.2f}  M={:.2f}".format(
                title_prefix, n_fish, burn_in_s + fi * dt, float(P_t[fi]), float(M_t[fi])
            )
        )
        return trails + bodies + [title]

    anim = animation.FuncAnimation(
        fig, update, frames=len(frame_ids), interval=int(1000 / fps), blit=False
    )
    anim.save(out_path, writer="pillow", fps=fps)
    plt.close(fig)
    return dict(P=float(P_t.mean()), M=float(M_t.mean()), path=out_path)


def make_video(
    out_path,
    phase="school",
    n_fish=100,
    k=1,
    gamma_att=None,
    gamma_ali=None,
    seed=71,
    dt=0.05,
    burn_in_s=100.0,
    record_s=25.0,
    fps=12,
):
    """Wang burst-and-coast S3/S4-style GIF."""
    if phase == "mill":
        gamma_att = 0.05 if gamma_att is None else gamma_att
        gamma_ali = 0.05 if gamma_ali is None else gamma_ali
    else:
        gamma_att = 0.04 if gamma_att is None else gamma_att
        gamma_ali = 0.20 if gamma_ali is None else gamma_ali

    t0 = _time.time()
    sim = WangBurstCoastSimulator(
        n=n_fish, k=k, gamma_att=gamma_att, gamma_ali=gamma_ali, seed=seed
    )
    print(
        "Wang {}: N={} k={} γ_Att={} γ_Ali={} burn={:.0f}s record={:.0f}s".format(
            phase, n_fish, k, gamma_att, gamma_ali, burn_in_s, record_s
        )
    )
    sim.run(dt, int(round(burn_in_s / dt)))
    n_rec = int(round(record_s / dt))
    xs = np.empty((n_rec, n_fish))
    ys = np.empty((n_rec, n_fish))
    phis = np.empty((n_rec, n_fish))
    for s in range(n_rec):
        sim.step(dt)
        xs[s], ys[s], phis[s] = sim.x, sim.y, sim.phi
    out = _render_gif(out_path, xs, ys, phis, dt, burn_in_s, fps, "Wang " + phase)
    print("  saved {} ({:.1f}s total)".format(out_path, _time.time() - t0))
    return out


def make_pid_video(
    out_path,
    n_fish=25,
    k=1,
    gamma_att=0.04,
    gamma_ali=0.20,
    seed=71,
    dt=0.05,
    burn_in_s=40.0,
    record_s=25.0,
    fps=12,
):
    """Stefan-like act()+PD tracking Wang setpoints."""
    t0 = _time.time()
    sim = StefanPIDSchool(
        n=n_fish, k=k, gamma_att=gamma_att, gamma_ali=gamma_ali, seed=seed
    )
    print(
        "StefanPID school: N={} k={} γ_Att={} γ_Ali={} burn={:.0f}s record={:.0f}s".format(
            n_fish, k, gamma_att, gamma_ali, burn_in_s, record_s
        )
    )
    sim.run(dt, int(round(burn_in_s / dt)))
    n_rec = int(round(record_s / dt))
    xs = np.empty((n_rec, n_fish))
    ys = np.empty((n_rec, n_fish))
    phis = np.empty((n_rec, n_fish))
    for s in range(n_rec):
        sim.step(dt)
        xs[s], ys[s], phis[s] = sim.x, sim.y, sim.phi
    out = _render_gif(out_path, xs, ys, phis, dt, burn_in_s, fps, "StefanPID")
    print("  saved {} ({:.1f}s total)".format(out_path, _time.time() - t0))
    return out


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    make_pid_video(
        os.path.join(here, "stefan_pid_school.gif"),
        n_fish=25,
        k=1,
        gamma_att=0.04,
        gamma_ali=0.20,
        seed=71,
        burn_in_s=40.0,
        record_s=25.0,
    )
