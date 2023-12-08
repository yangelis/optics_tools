# %%
from misc import *

try:
    get_ipython().run_line_magic("matplotlib", "inline")
except:
    pass


# %%

hpre = os.getenv("HOME")
prefix = f"{hpre}/Projects/hllhc_optics/"

f1 = prefix + "summer_studies/collapse/opt_collapse_flathv_900_1800_1500_thin.madx"
simple_name = f1.split("/")[-1].split(".")[0] + ".json"


collider = None
if os.path.exists(prefix + simple_name):
    collider = xt.Multiline.from_json(prefix + simple_name)
else:
    collider = build_collider_from_mad(f1, prefix=prefix)
    collider.to_json(prefix + simple_name)

# %%
collider.lhcb1.particle_ref = xp.Particles(p0c=7000e9, q0=1, mass0=xp.PROTON_MASS_EV)
collider.lhcb2.particle_ref = xp.Particles(p0c=7000e9, q0=1, mass0=xp.PROTON_MASS_EV)
collider.build_trackers()

# %%
collider.lhcb1.twiss_default["method"] = "4d"
collider.lhcb2.twiss_default["method"] = "4d"
collider.lhcb2.twiss_default["reverse"] = True
# %%
my_twiss = collider.twiss(method="4d")
print(my_twiss["lhcb1"][["x", "y", "betx", "bety", "px", "py"], ["ip1", "ip5"]])
print(my_twiss["lhcb2"][["x", "y", "betx", "bety", "px", "py"], ["ip1", "ip5"]])


# %%
collider["lhcb1"].varval["i_oct_b1"] = 450
collider["lhcb1"].varval["on_x1"] = 250
collider["lhcb1"].varval["on_x5"] = 250
line = collider["lhcb1"].cycle("ip3")
# %%


zeta0 = 0
delta0 = 0

n_x_norm = 10
n_y_norm = 10
# n_x_norm = 2
# n_y_norm = 1

nemitt_x = 2.3e-6
nemitt_y = 2.3e-6
num_turns = 100
zeropad = 100000

x_norm_range = (0.1, 6)
y_norm_range = (0.1, 6)

Jx_min = nemitt_x * x_norm_range[0] ** 2 / 2
Jx_max = nemitt_x * x_norm_range[1] ** 2 / 2
Jy_min = nemitt_y * y_norm_range[0] ** 2 / 2
Jy_max = nemitt_y * y_norm_range[1] ** 2 / 2

Jx_grid = np.linspace(Jx_min, Jx_max, n_x_norm)
Jy_grid = np.linspace(Jy_min, Jy_max, n_y_norm)

Jx_2d, Jy_2d = np.meshgrid(Jx_grid, Jy_grid)

x_norm_2d = np.sqrt(2 * Jx_2d / nemitt_x)
y_norm_2d = np.sqrt(2 * Jy_2d / nemitt_y)

# %%

fig, axs = plt.subplots(figsize=(8, 8))
axs.plot(Jx_2d, Jy_2d, "o", color="b")
axs.set_xlabel(r"$J_x$")
axs.set_ylabel(r"$J_y$")
plt.show()

# %%

particles = line.build_particles(
    x_norm=x_norm_2d.flatten(),
    y_norm=y_norm_2d.flatten(),
    nemitt_x=nemitt_x,
    nemitt_y=nemitt_y,
    zeta=zeta0,
    delta=delta0,
    freeze_longitudinal=True,
    method="4d",
)

# %%
line.track(
    particles, num_turns=num_turns, turn_by_turn_monitor=True, freeze_longitudinal=True
)

x = line.record_last_track.x
y = line.record_last_track.y

# %%
ipar = 0
fig, axs = plt.subplots(figsize=(21, 11))
axs.plot(x[ipar], y[ipar], "o", color="black")
axs.set_xlabel("x")
axs.set_ylabel("y")
axs.set_title(f"Particle {ipar}")
plt.show()
# %%

fig, axs = plt.subplots(figsize=(21, 11))
for ipar in range(x.shape[0]):
    axs.plot(x[ipar], y[ipar], "o", label=f"Particle {ipar}")
axs.set_xlabel("x")
axs.set_ylabel("y")
# axs.legend()
plt.show()

# %%
fig, axs = plt.subplots(figsize=(21, 11))
axs.plot(x, y, "o", color="black")
axs.set_xlabel("x")
axs.set_ylabel("y")
plt.show()
# %%
mon = line.record_last_track
# %%
x_noCO = mon.x - np.atleast_2d(mon.x.mean(axis=1)).T
y_noCO = mon.y - np.atleast_2d(mon.y.mean(axis=1)).T

# %%
fig, axs = plt.subplots(2, 1, figsize=(21, 11))
axs[0].plot(mon.x)
axs[0].set_xlabel("N Particle")
axs[0].set_ylabel("x")
axs[1].plot(mon.y)
axs[1].set_xlabel("N Particle")
axs[1].set_ylabel("y")
plt.show()

# %%
fig, axs = plt.subplots(2, 1, figsize=(21, 11))
axs[0].plot(mon.x.T)
# axs[0].plot(x_noCO.T)
axs[0].set_xlabel("N Turns")
# axs[0].set_ylabel(r"$x_{co}$")
axs[0].set_ylabel("x")
axs[1].plot(mon.y.T)
# axs[1].plot(y_noCO.T)
axs[1].set_xlabel("N Turns")
# axs[1].set_ylabel(r"$y_{co}$")
axs[1].set_ylabel("y")
plt.show()

# %%
frequency = np.fft.fftfreq(zeropad)
qx = [
    np.abs(
        frequency[
            np.argmax(np.abs(np.fft.fft(x[ii] * np.hanning(x.shape[1]), n=zeropad)))
        ]
    )
    for ii in range(x.shape[0])
]
qy = [
    np.abs(
        frequency[
            np.argmax(np.abs(np.fft.fft(y[ii] * np.hanning(x.shape[1]), n=zeropad)))
        ]
    )
    for ii in range(x.shape[0])
]

qx = np.reshape(qx, x_norm_2d.shape)
qy = np.reshape(qy, y_norm_2d.shape)

# %%
fig, axs = plt.subplots(figsize=(21, 11))
axs.plot(qx, qy, "o-", color="black")
axs.plot(qx.T, qy.T, "o-", color="black")
axs.set_xlabel(r"$q_x$")
axs.set_ylabel(r"$q_x$")
plt.show()

# %%
### Analytical

tw1 = collider.lhcb1.twiss(method="4d", strengths=True)


oct_tw = tw1[:, "mo.*"]

dqxdex = 1 / 32 / np.pi * np.sum(oct_tw.betx**2 * oct_tw.k3nl)

dqxdey = -1 / 16 / np.pi * np.sum(oct_tw.betx * oct_tw.bety * oct_tw.k3nl)

dqydey = 1 / 32 / np.pi * np.sum(oct_tw.bety**2 * oct_tw.k3nl)

dqydex = dqxdey

tw_quads = tw1[:, "mq.*"]

quad_dqx = 1 / 4 / np.pi * np.sum(tw_quads.betx * tw_quads.k1nl)
quad_dqy = -1 / 4 / np.pi * np.sum(tw_quads.bety * tw_quads.k1nl)


# %%
qx_an = tw1.qx + dqxdex * (2 * Jx_2d) + dqxdey * (2 * Jy_2d)
qy_an = tw1.qy + dqydey * (2 * Jy_2d) + dqydex * (2 * Jx_2d)


jx = Jx_grid / collider["lhcb1"].particle_ref.gamma0
jy = Jy_grid / collider["lhcb1"].particle_ref.gamma0

xv, yv = np.meshgrid(jx, jy)
positions = np.vstack([xv.ravel(), yv.ravel()]).T

qx_an = tw1.qx + dqxdex * (2 * xv) + dqxdey * (2 * yv)  # + quad_dqx
qy_an = tw1.qy + dqydey * (2 * yv) + dqydex * (2 * xv)  # + quad_dqy


# %%
fig, axs = plt.subplots(figsize=(21, 11))
axs.plot(qx.T, qy.T, "o-", color="blue")
axs.plot(qx, qy, "o-", color="blue")
axs.plot(qx.flatten(), qy.flatten(), "o-", color="blue", label="Tracking")

axs.plot(qx_an - 62.0, qy_an - 60.0, "o-", color="black")
axs.plot((qx_an - 62.0).T, (qy_an - 60.0).T, "o-", color="black")
axs.plot(
    (qx_an - 62.0).flatten(),
    (qy_an - 60.0).flatten(),
    "o-",
    color="black",
    label="Analytical",
)

axs.set_xlabel(r"$q_x$")
axs.set_ylabel(r"$q_y$")
axs.legend()
plt.show()
# %%
# normalized emmitance
egeom_x = nemitt_x / line.particle_ref.gamma0
egeom_y = nemitt_y / line.particle_ref.gamma0


fig, axs = plt.subplots(2, 1, figsize=(21, 11))
# axs[0].plot(xv, qx)
# axs[1].plot(yv, qy)
axs[0].plot(xv.flatten(), qx.flatten(), "o")
axs[1].plot(yv.flatten(), qy.flatten(), "o")
# axs[0].plot(Jx_2d * egeom_x, qx)
# axs[1].plot(Jy_2d * egeom_y, qy)
# axs[0].plot(Jx_2d.flatten() * egeom_x, qx.flatten())
# axs[1].plot(Jy_2d.flatten() * egeom_y, qy.flatten())

plt.show()

# %%
fxx = np.polyfit(xv.flatten(), qx.flatten(), 1)
# %%
z1 = np.poly1d(fxx)
dqxdex_track = z1[1]
# %%
fig, axs = plt.subplots(2, 1, figsize=(21, 11))

axs[0].plot(xv.flatten(), qx.flatten(), "o")
axs[0].plot(xv.flatten(), z1(xv.flatten()))

plt.show()

# %%


fig, axs = plt.subplots(figsize=(21, 11))
axs.plot(xv.flatten(), qx.flatten(), "o")
# axs.plot(yv.flatten(), qy.flatten(), "o")
plt.show()
fxx = np.polyfit(xv.flatten(), qx.flatten(), 1)
z1 = np.poly1d(fxx)
dqxdex_track = z1[1]
fig, axs = plt.subplots(figsize=(21, 11))
axs.plot(xv.flatten(), qx.flatten(), "o")
axs.plot(xv.flatten(), z1(xv.flatten()))
plt.show()
print(dqxdex_track / 2, dqxdex)

# %%
dqxdex_tracks = []
fits = []
for i in range(xv.shape[0]):
    dqxdex_tracks.append(np.polyfit(xv[i, :], qx[i, :], 1)[0])
    fits.append(np.polyfit(xv[i, :], qx[i, :], 1))

# %%
fig, axs = plt.subplots(figsize=(21, 11))
for i in range(xv.shape[0]):
    fff = np.poly1d(fits[i])
    axs.plot(xv.flatten(), fff(xv.flatten()))

plt.show()


# %%
print(np.asarray(dqxdex_tracks) / 2)
mean_dqxdex_tracks = np.mean(dqxdex_tracks) / 2
print(mean_dqxdex_tracks)
# %%
