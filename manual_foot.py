# %%
from misc import *

try:
    get_ipython().run_line_magic("matplotlib", "inline")
except:
    pass


# %%


# %%

hpre = os.getenv("HOME")
prefix = f"{hpre}/Projects/hllhc_optics/"

f1 = (
    prefix
    # + "opt_flathv_500_2000_thin_testing.madx"
    # + "acc-models-lhc/strengths/round/levelling/opt_levelling_580_1500_thin.madx"
    # + "acc-models-lhc/strengths/flat/opt_flathv_500_2000_thin.madx"
    # + "summer_studies/collapse/opt_collapse_700_1500_thin.madx"
    # + "summer_studies/collapse/opt_collapse_1100_1500_thin.madx"
    # + "summer_studies/collapse/opt_collapse_1000_1500_no_triplet_match.madx",
    # + "summer_studies/collapse/opt_collapse_1000_1500_rematch_thin.madx"
    # + "summer_studies/collapse/opt_collapse_1000_1500_no_triplet_match_thin.madx"
    # + "summer_studies/collapse/opt_collapse_1000_1500_triplet_match3_thin_again.madx"
    # + "summer_studies/collapse/opt_flathv_450_1800_1500_thin.madx"
    # + "start_of_collapse/flat/opt_collapse_flathv_700_2800_test.madx"
    # + "summer_studies/collapse/opt_collapse_700_1500_rematched_again_thin.madx"
    # + "summer_studies/collapse/opt_collapse_flathv_700_2800_thin.madx"
    + "summer_studies/collapse/opt_collapse_flathv_900_1800_1500_thin.madx"
    # + "summer_studies/flat_temp/reduce/opt_flat_606_1212_1500_step.madx"
    # + "summer_studies/collapse/opt_flathv_450_1800_1500.madx"
    # + "summer_studies/flat_temp/opt_flathv_455_910_1500.madx"
)
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
print(my_twiss.lhcb1[["x", "y", "betx", "bety", "px", "py"], ["ip1", "ip5"]])
print(my_twiss.lhcb2[["x", "y", "betx", "bety", "px", "py"], ["ip1", "ip5"]])


# %%
line = collider["lhcb1"]

nemitt_x = 2.3e-6
nemitt_y = 2.3e-6
num_turns = 1000
zeropad = 100000

# %%
"""
returns axx, axy, ayx, ayy that correspond to the first order amplitude-detuning,
estimated through tracking.
The relations between tune (Q) and action (J) in this case is:
$$ Q_x = Q_{x,0} + axx * Jx + axy * Jy $$
$$ Q_y = Q_{y,0} + ayx * Jx + ayy * Jy $$
"""

# normalized emmitance
egeom_x = nemitt_x / line.particle_ref.gamma0
egeom_y = nemitt_y / line.particle_ref.gamma0

frequency = np.fft.fftfreq(zeropad)

# %%

sigma = 6
num_r = 100
JJ = np.linspace(0.01, sigma, num_r)
A_norm = np.sqrt(2 * JJ).flatten()
other_norm = 0.01

particles = line.build_particles(
    x_norm=A_norm,
    y_norm=other_norm,
    nemitt_x=nemitt_x,
    nemitt_y=nemitt_y,
    num_particles=num_r,
)
# %%
line.track(particles, num_turns=num_turns, turn_by_turn_monitor=True)
# line.tracker._context.synchronize()
# %%
x = line.record_last_track.x
y = line.record_last_track.y
# %%
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

# %%
fxx = np.polyfit(JJ * egeom_x, qx, 1)
fyx = np.polyfit(JJ * egeom_x, qy, 1)
axx = np.polyfit(JJ * egeom_x, qx, 1)[0]
ayx = np.polyfit(JJ * egeom_x, qy, 1)[0]

z1 = np.poly1d(fxx)
z2 = np.poly1d(fyx)


fig, axs = plt.subplots(2, 1, figsize=(21, 11))
axs[0].plot(JJ * egeom_x, qx)
axs[0].plot(JJ * egeom_x, z1(JJ * egeom_x))
axs[1].plot(JJ * egeom_x, qy)
axs[1].plot(JJ * egeom_x, z2(JJ * egeom_x))
plt.show()


# %%
# switch x and y
particles = line.build_particles(
    x_norm=other_norm,
    y_norm=A_norm,
    nemitt_x=nemitt_x,
    nemitt_y=nemitt_y,
    num_particles=num_r,
)

line.track(particles, num_turns=num_turns, turn_by_turn_monitor=True)
# line.tracker._context.synchronize()

x = line.record_last_track.x
y = line.record_last_track.y

qx = [
    abs(
        frequency[
            np.argmax(np.abs(np.fft.fft(x[ii] * np.hanning(x.shape[1]), n=zeropad)))
        ]
    )
    for ii in range(x.shape[0])
]
qy = [
    abs(
        frequency[
            np.argmax(np.abs(np.fft.fft(y[ii] * np.hanning(x.shape[1]), n=zeropad)))
        ]
    )
    for ii in range(x.shape[0])
]

axy = np.polyfit(JJ * egeom_y, qx, 1)[0]
ayy = np.polyfit(JJ * egeom_y, qy, 1)[0]

fxy = np.polyfit(JJ * egeom_y, qx, 1)
fyy = np.polyfit(JJ * egeom_y, qy, 1)
# %%
z3 = np.poly1d(fxy)
z4 = np.poly1d(fyy)


fig, axs = plt.subplots(2, 1, figsize=(21, 11))
axs[0].plot(JJ * egeom_y, qx)
axs[0].plot(JJ * egeom_y, z3(JJ * egeom_y))
axs[1].plot(JJ * egeom_y, qy)
axs[1].plot(JJ * egeom_y, z4(JJ * egeom_y))
plt.show()
# %%


tw1 = collider["lhcb1"].twiss(method="4d", strengths=True)
array_Jx = np.linspace(0.01, 6, 100)
array_Jy = np.linspace(0.01, 6, 100)


qx_t = tw1.qx + axx * (2 * array_Jx) + axy * (2 * array_Jy)
qy_t = tw1.qy + ayy * (2 * array_Jy) + ayx * (2 * array_Jx)
# %%
# print(qx_t, qy_t)

fig, axs = plt.subplots()
axs.plot(qx_t, qy_t, "o")
plt.show()
# %%

qx_t = tw1.qx + axx * (array_Jx) + axy * (array_Jy)
qy_t = tw1.qy + ayy * (array_Jy) + ayx * (array_Jx)
print(qx_t, qy_t)
fig, axs = plt.subplots()
axs.plot(qx_t, qy_t, "o")
plt.show()

# %%
