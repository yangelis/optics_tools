# %%
from misc import *

try:
    get_ipython().run_line_magic("matplotlib", "inline")
except:
    pass

# %%
hpre = os.getenv("HOME")
prefix = f"{hpre}/Projects/hllhc_optics/"

optics_files = [
    prefix + "summer_studies/collapse/opt_collapse_flathv_900_1800_1500_thin.madx",
    # prefix + "summer_studies/collapse/opt_collapse_flathv_700_2800_thin.madx",
    # prefix + "summer_studies/collapse/opt_flathv_600_1200_1500_thin.madx",
    # prefix + "summer_studies/collapse/opt_flathv_450_1800_1500_thin.madx",
]

# %%

results = []
colliders = []

for optics_file in optics_files:
    simple_name = optics_file.split("/")[-1].split(".")[0] + ".json"

    collider = None
    if os.path.exists(prefix + simple_name):
        collider = xt.Multiline.from_json(prefix + simple_name)
    else:
        collider = build_collider_from_mad(optics_file, prefix=prefix)
        collider.to_json(prefix + simple_name)

    collider.lhcb1.twiss_default["method"] = "4d"
    collider.lhcb2.twiss_default["method"] = "4d"
    collider.lhcb2.twiss_default["reverse"] = True

    collider.lhcb1.particle_ref = xp.Particles(
        p0c=7000e9, q0=1, mass0=xp.PROTON_MASS_EV
    )
    collider.lhcb2.particle_ref = xp.Particles(
        p0c=7000e9, q0=1, mass0=xp.PROTON_MASS_EV
    )
    collider.build_trackers()

    colliders.append(collider)


# %%
# mos = np.round(np.linspace(-400, 0, 500), decimals=2)
mos = np.round(np.arange(-400, 500, 50), decimals=2)

for col in colliders:
    dets = detuning_terms(col, mo_list=mos)

    results.append(dets)


# %%
def detuning_terms3(collider, mo_list, bim=1):
    results = {
        "dqxdex": [],
        "dqxdey": [],
        "dqydey": [],
    }

    for mo in mo_list:
        # collider[f"lhcb{bim}"].varval[f"i_oct_b{bim}"] = mo

        for ia in ["23", "34", "67", "78"]:
            collider.vars[f"kof.a{ia}b1"] = 60
            collider.vars[f"kod.a{ia}b1"] = 60
        for ia in ["81", "12", "45", "56"]:
            collider.vars[f"kof.a{ia}b1"] = mo
            collider.vars[f"kod.a{ia}b1"] = mo

        tw1 = collider[f"lhcb{bim}"].twiss(method="4d", strengths=True)

        arc_octupoles = []
        for ii, (i, j) in enumerate(arcs):
            oct_arc = tw1[:, f"s.ds.r{i}.b{bim}":f"e.ds.l{j}.b{bim}"][:, "mo.*"]
            arc_octupoles.append(oct_arc)

        oct_tw = xt.TwissTable.concatenate(arc_octupoles)

        dqxdex = 1 / 32 / np.pi * np.sum(oct_tw.betx**2 * oct_tw.k3nl)

        dqxdey = -1 / 16 / np.pi * np.sum(oct_tw.betx * oct_tw.bety * oct_tw.k3nl)

        dqydey = 1 / 32 / np.pi * np.sum(oct_tw.bety**2 * oct_tw.k3nl)

        results["dqxdex"].append(dqxdex)
        results["dqxdey"].append(dqxdey)
        results["dqydey"].append(dqydey)
    collider[f"lhcb{bim}"].varval[f"i_oct_b{bim}"] = 0
    return results


# %%
mos = np.round(np.arange(-400, 500, 50), decimals=2)
results3 = []
for col in colliders:
    dets = detuning_terms3(col, mo_list=mos)

    results3.append(dets)

# %%

# for ia in lm.ARC_NAMES:
#     # colliders[0].varval["i_oct_b1"] = 100
#     print(colliders[0].vars[f"kof.a{ia}b1"]._expr())
#     print(colliders[0].vars[f"kod.a{ia}b1"]._expr())
#    # colliders[0].varval["i_oct_b1"] = 0
# %%
labels = [
    r"$\beta_x=1.8 [m],\beta_y=0.9 [m]$",
    r"$\beta_x=2.8 [m],\beta_y=0.7 [m]$",
    r"$\beta_x=1.2 [m],\beta_y=0.6 [m]$",
    r"$\beta_x=1.8 [m],\beta_y=0.45[m]$",
]

# %%
fig, axs = plt.subplots(1, 3, figsize=(31, 11))
for i, res in enumerate(results):
    axs[0].plot(mos, np.array(res["dqxdex"]) * 1e-5, "o-", label=labels[i])

    axs[1].plot(mos, np.array(res["dqydey"]) * 1e-5, "o-", label=labels[i])

    axs[2].plot(mos, np.array(res["dqxdey"]) * 1e-5, "o-", label=labels[i])


for i in range(3):
    axs[i].grid(True)
    axs[i].set_xlabel("Octupole Current [A]")
    axs[i].legend()

axs[0].set_ylabel(r"$\partial Q_x/\partial \epsilon_x [m^{-1}] 10^5$")
axs[1].set_ylabel(r"$\partial Q_y /\partial \epsilon_y [m^{-1}] 10^5$")
axs[2].set_ylabel(r"$\partial Q_x /\partial \epsilon_y [m^{-1}] 10^5$")


plt.show()
# %%
fig, axs = plt.subplots(1, 3, figsize=(31, 11))
for i, (res, res3) in enumerate(zip(results, results3)):
    axs[0].plot(mos, np.array(res["dqxdex"]) * 1e-5, "o-", label=labels[i])
    axs[0].plot(
        mos, np.array(res3["dqxdex"]) * 1e-5, "o-", label=labels[i], color="red"
    )

    axs[1].plot(mos, np.array(res["dqydey"]) * 1e-5, "o-", label=labels[i])
    axs[1].plot(
        mos, np.array(res3["dqydey"]) * 1e-5, "o-", label=labels[i], color="red"
    )

    axs[2].plot(mos, np.array(res["dqxdey"]) * 1e-5, "o-", label=labels[i])
    axs[2].plot(
        mos, np.array(res3["dqxdey"]) * 1e-5, "o-", label=labels[i], color="red"
    )


for i in range(3):
    axs[i].grid(True)
    axs[i].set_xlabel("Octupole Current [A]")
    axs[i].legend()

axs[0].set_ylabel(r"$\partial Q_x/\partial \epsilon_x [m^{-1}] 10^5$")
axs[1].set_ylabel(r"$\partial Q_y /\partial \epsilon_y [m^{-1}] 10^5$")
axs[2].set_ylabel(r"$\partial Q_x /\partial \epsilon_y [m^{-1}] 10^5$")


plt.show()

# %%

iocts = np.arange(0, 600, 100)
# ic = 3

# fps = []

# temp_tw = colliders[ic]["lhcb1"].twiss(method="4d")[["betx", "bety"], "ip1"]
# betx, bety = temp_tw.betx[0], temp_tw.bety[0]

# for i in range(iocts.size):
#     colliders[ic]["lhcb1"].varval["i_oct_b1"] = iocts[i]
#     lb1 = colliders[ic]["lhcb1"].cycle("ip3")
#     fp0 = lb1.get_footprint(
#         nemitt_x=2.3e-6,
#         nemitt_y=2.3e-6,
#         freeze_longitudinal=True,
#         mode="polar",
#     )
#     fps.append(fp0)

# # %%

cmap = plt.get_cmap("inferno_r")
colors = [cmap(i) for i in np.linspace(0.1, 0.8, iocts.size)]

# fig, axs = plt.subplots(figsize=(11, 11))
# for i, fp0 in enumerate(fps):
#     fp0.plot(
#         color=colors[i],
#         label=r"$I_{oct}" + f"={iocts[i]} A$",
#         ax=axs,
#         marker="o",
#         ms=4,
#         linestyle="-.",
#         linewidth=0.5,
#     )
# axs.legend()
# axs.set_title(rf"$\beta_x={betx:1.2f}\ m, \beta_y={bety:1.2f}\ m$")

# plt.show()


# %%
def footprints(colliders, ioct):
    cfps = []
    bets = []

    for ic in range(len(colliders)):
        colliders[ic]["lhcb1"].varval["i_oct_b1"] = ioct
        lb1 = colliders[ic]["lhcb1"].cycle("ip3")
        fp0 = lb1.get_footprint(
            nemitt_x=2.3e-6,
            nemitt_y=2.3e-6,
            freeze_longitudinal=True,
            mode="polar",
        )
        temp_tw = colliders[ic]["lhcb1"].twiss(method="4d")[["betx", "bety"], "ip1"]
        betx, bety = temp_tw.betx[0], temp_tw.bety[0]
        cfps.append(fp0)
        bets.append((betx, bety))

    fig, axs = plt.subplots(figsize=(11, 11))
    for i, fp0 in enumerate(cfps):
        betx, bety = bets[i]
        fp0.plot(
            color=colors[i],
            label=rf"$\beta_x={betx:1.2f}\ m, \beta_y={bety:1.2f}\ m$",
            ax=axs,
            marker="o",
            ms=4,
            linestyle="-.",
            linewidth=0.5,
        )
    axs.legend()
    axs.set_title(
        rf"$I_{{oct}}={ioct} A$",
    )

    plt.show()


# %%
# for mo in mos:
# footprints(colliders, mo)

# %%


def detuning_terms2(collider, mo, jx, jy, bim=1):
    results = {"dqxdex": None, "dqxdey": None, "dqydey": None, "dqx": None, "qx": None}

    collider[f"lhcb{bim}"].varval[f"i_oct_b{bim}"] = mo

    tw1 = collider[f"lhcb{bim}"].twiss(method="4d", strengths=True)

    arc_octupoles = []
    for ii, (i, j) in enumerate(arcs):
        oct_arc = tw1[:, f"s.ds.r{i}.b{bim}":f"e.ds.l{j}.b{bim}"][:, "mo.*"]
        arc_octupoles.append(oct_arc)

    oct_tw = xt.TwissTable.concatenate(arc_octupoles)

    dqxdex = 1 / 32 / np.pi * np.sum(oct_tw.betx**2 * oct_tw.k3nl)

    dqxdey = -1 / 16 / np.pi * np.sum(oct_tw.betx * oct_tw.bety * oct_tw.k3nl)

    dqydey = 1 / 32 / np.pi * np.sum(oct_tw.bety**2 * oct_tw.k3nl)

    dqydex = dqxdey

    results["dqxdex"] = dqxdex
    results["dqxdey"] = dqxdey
    results["dqydey"] = dqydey
    results["qx_rms"] = dqxdex**2 * (2 * jx) ** 2 + dqxdey**2 * (2 * jy) ** 2
    results["dqx"] = dqxdex * (2 * jx) + dqxdey * (2 * jy)
    results["qx"] = tw1.qx + dqxdex * (2 * jx) + dqxdey * (2 * jy)
    results["qy"] = tw1.qy + dqydey * (2 * jy) + dqydex * (2 * jx)
    collider[f"lhcb{bim}"].varval[f"i_oct_b{bim}"] = 0
    return results


# %%
detuning_terms2(colliders[0], 450, 0, 0)

# %%
jx = 6 * np.linspace(0.01, 6, 100) * 1e-6 / colliders[0]["lhcb1"].particle_ref.gamma0
jy = 6 * np.linspace(0.01, 6, 100) * 1e-6 / colliders[0]["lhcb1"].particle_ref.gamma0

xv, yv = np.meshgrid(jx, jy)
positions = np.vstack([xv.ravel(), yv.ravel()]).T


# %%
jjx = positions[:, 0]
jjy = positions[:, 1]
# plt.plot(jjx, jjy, "p")
res = detuning_terms2(colliders[0], 450, jjx, jjy)
# %%
fig, axes = plt.subplots(figsize=(21, 8))
axes.plot(res["qx"], res["qy"], "o")
plt.show()
# %%
jjx = np.linspace(0, 1, 10) * 1e-6
rr = np.zeros(shape=jjx.size)
res2 = []
for i, ix in enumerate(jjx):
    res2.append(detuning_terms2(colliders[0], 450, ix, 2.5e-6 / 7463))
    rr[i] = res2[-1]["dqx"]

print(rr)

# %%
qxs = np.fromiter(map(lambda d: d["qx"], res2), dtype=np.float64)

# %%
plt.plot(jjx, qxs, marker="o", linestyle="none")
plt.show()

# %%
rr = np.zeros(shape=positions.shape[0])
res2 = []
for i, (ix, iy) in enumerate(positions):
    res2.append(detuning_terms2(colliders[0], 450, ix, iy))


# %%
qxs2 = np.fromiter(map(lambda d: d["qx"], res2), dtype=np.float64)
qys2 = np.fromiter(map(lambda d: d["qy"], res2), dtype=np.float64)

# %%
plt.plot(positions[:, 0], np.sqrt(rr), marker="o", linestyle="none")
plt.show()
# %%

# %%
hw, yedges, xedges = np.histogram2d(
    xv.flatten(), yv.flatten(), weights=qxs2 - 62.30, bins=(5, 5)
)
fig, axes = plt.subplots(figsize=(21, 8))
pcm = axes.pcolormesh(xedges, yedges, hw.T, rasterized=True)
fig.colorbar(pcm, ax=axes, pad=0)
plt.show()

# %%
plt.plot(qys2 - 60.32, marker="o", linestyle="")
plt.show()
# %%
tw0 = colliders[0].lhcb1.twiss(method="4d")[:, "mo.*"]
# %%
fig, axs = plt.subplots(2, 1, figsize=(21, 20))
for col in colliders:
    tw_temp = col["lhcb1"].twiss(method="4d")
    tw_temp_mo = tw_temp[:, "mo.*"]
    axs[0].plot(
        tw_temp_mo.s,
        tw_temp_mo.betx,
        label=f"{np.min(tw_temp.betx):2.2f}",
        marker="o",
        linestyle="",
    )
    axs[1].plot(
        tw_temp_mo.s,
        tw_temp_mo.bety,
        label=f"{np.min(tw_temp.betx):2.2f}",
        marker="o",
        linestyle="",
    )
axs[0].set_xlim(0, 5000)
axs[1].set_xlim(0, 5000)

axs[0].legend()
axs[1].legend()
plt.show()

# %%

colliders[0]["lhcb1"].varval["i_oct_b1"] = 450
lb1 = colliders[0]["lhcb1"].cycle("ip3")
fp0 = lb1.get_footprint(
    nemitt_x=2.3e-6,
    nemitt_y=2.3e-6,
    freeze_longitudinal=True,
    delta0=0,
    zeta0=0,
    mode="uniform_action_grid",
)

# %%
fig, axs = plt.subplots(figsize=(21, 11))
fp0.plot(
    ax=axs,
    marker="o",
    ms=4,
    linestyle="-.",
    linewidth=0.5,
)
axs.plot(res["qx"] - 62.0, res["qy"] - 60.0, "p")
plt.show()

# %%

# %%
fig, axs = plt.subplots(figsize=(18, 18))
axs.set_title("Beam 1")
axs.plot(tw_k["lhcb1"]["s"], tw_k["lhcb1"]["x"], label="x", color="black", lw=3)
axs.plot(tw_k["lhcb1"]["s"], tw_k["lhcb1"]["y"], label="y", color="red", lw=3)
axs00 = axs.twinx()
axs00.plot(tw_k["lhcb1"]["s"], tw_k["lhcb1"]["dx"], label="dx")
axs00.plot(tw_k["lhcb1"]["s"], tw_k["lhcb1"]["dy"], label="dy")
axs00.plot(
    tw_k["lhcb1"][["s"], "ms.*"],
    [2] * tw_k["lhcb1"][["s"], "ms.*"].size,
    "o",
    ms=6,
    label="Octupoles",
)
axs00.plot(
    tw_k["lhcb1"][["s"], "mo.*"],
    [2] * tw_k["lhcb1"][["s"], "mo.*"].size,
    "o",
    ms=6,
    label="Octupoles",
)
axs.legend(loc=2)
axs00.legend(loc=1)
axs.set_xlabel("s [m]")
axs.set_ylabel(r"$co [m]$")
axs00.set_ylabel("Dx,Dy [m]")
plt.show()


# %%
fig, axs = plt.subplots(figsize=(10, 10))
iocts = [450, 300, 200, 150, -150, -200, -300, -450]
cmap = plt.get_cmap("viridis")
colors = [cmap(i) for i in np.linspace(0, 1, len(iocts))]

fps = []
for i in range(4):
    lb1.vars["i_oct_b1"] = iocts[i]
    fp0 = lb1.get_footprint(
        nemitt_x=2.3e-6,
        nemitt_y=2.3e-6,
        freeze_longitudinal=True,
        mode="polar",
    )
    fp0.plot(
        color=colors[i],
        label=f"I_oct={iocts[i]} A",
        ax=axs,
        marker="o",
        ms=3,
        linestyle="--",
    )
    fps.append(fp0)
axs.legend()

plt.show()

# %%
collider["lhcb1"].varval["i_oct_b1"] = 0
# %%


collider["lhcb1"].varval["i_oct_b1"] = 450
tw1 = collider["lhcb1"].twiss(method="4d", strengths=True)


# %%

all_oct = tw_k["lhcb1"][:, "mo.*"]
fig, axs = plt.subplots(figsize=(11, 11))
axs.plot(all_oct.s, all_oct.betx, color="blue")
axs.plot(all_oct.s, all_oct.bety, color="red")
plt.show()

# %%
arc_octupoles = []

for ii, (i, j) in enumerate(arcs):
    oct_arc = tw_k["lhcb1"][:, f"s.ds.r{i}.b{bim}":f"e.ds.l{j}.b{bim}"][:, "mo.*"]
    arc_octupoles.append(oct_arc)

print(arc_octupoles[0]["k3nl"])
print(arc_octupoles[0].name)

oct_tw = xt.TwissTable.concatenate(arc_octupoles)

# %%

s1 = np.sum(
    [
        arc_octupoles[i].betx ** 2 * arc_octupoles[i].k3nl
        for i in range(len(arc_octupoles))
    ]
)

dqxdex = 1 / 32 / np.pi * np.sum(oct_tw.betx**2 * oct_tw.k3nl)
# dqxdex = 1 / 32 / np.pi * float(s1)

dqxdey = -1 / 16 / np.pi * np.sum(oct_tw.betx * oct_tw.bety * oct_tw.k3nl)

dqydey = 1 / 32 / np.pi * np.sum(oct_tw.bety**2 * oct_tw.k3nl)

print(f"{dqxdex*1e-5=}, {dqxdey*1e-5=}, {dqydey*1e-5=}")


# %%
print(
    f"{np.sqrt((dqxdex * 2.5e-6/7463) ** 2 + (dqxdey * 2.5e-6/7463) ** 2)*1e4} * 1e-4"
)

# %%
fig, axs = plt.subplots(figsize=(11, 11))
for i in range(len(arc_octupoles)):
    axs.plot(arc_octupoles[i].s, arc_octupoles[i].betx, color="blue")
    axs.plot(arc_octupoles[i].s, arc_octupoles[i].bety, color="red")
plt.show()


# %%


# %%
mos = np.round(np.linspace(-400, 0, 20), decimals=2)
a = detuning_terms(collider, mos)

# %%
fig, axs = plt.subplots(1, 3, figsize=(31, 11))
axs[0].plot(mos, np.array(a["dqxdex"]) * 1e-5, "o-")

axs[1].plot(mos, np.array(a["dqydey"]) * 1e-5, "o-")

axs[2].plot(mos, np.array(a["dqxdey"]) * 1e-5, "o-")


axs[0].grid(True)
axs[1].grid(True)
axs[2].grid(True)

axs[0].set_xlabel("Octupole Current [A]")
axs[1].set_xlabel("Octupole Current [A]")
axs[2].set_xlabel("Octupole Current [A]")

axs[0].set_ylabel(r"$\partial Q_x/\partial \epsilon_x [m^{-1}] 10^5$")
axs[1].set_ylabel(r"$\partial Q_y /\partial \epsilon_y [m^{-1}] 10^5$")
axs[2].set_ylabel(r"$\partial Q_x /\partial \epsilon_y [m^{-1}] 10^5$")


plt.show()
# %%
