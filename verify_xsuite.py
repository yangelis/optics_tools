# %%
from misc import *

try:
    get_ipython().run_line_magic("matplotlib", "inline")
except:
    pass


# %%
import ruamel.yaml

ryaml = ruamel.yaml.YAML()


# %%
def plot_arcs_xy(collider, bim=1):
    fig, axs = plt.subplots(8, figsize=(10, 60))
    for ii, (i, j) in enumerate(arcs):
        tw_arcij = twiss_arc(collider, i, j, bim)
        axs[ii].set_title(f"Arc {i}{j}")
        axs[ii].plot(tw_arcij["s"], tw_arcij["x"], label="x")
        axs[ii].plot(tw_arcij["s"], tw_arcij["y"], label="y")
        axs[ii].set_xlabel("s [m]")
        axs[ii].set_ylabel(r"$co [m]$")
        axs[ii].legend()
        # axs[ii].set_ylim(-0.000001, 0.000001)
    plt.show()


def plot_arcs_betas(collider, bim=1):
    fig, axs = plt.subplots(8, figsize=(10, 60))
    for ii, (i, j) in enumerate(arcs):
        tw_arcij = twiss_arc(collider, i, j, bim)
        axs[ii].set_title(f"Arc {i}{j}")
        axs[ii].plot(tw_arcij["s"], tw_arcij["betx"], label=r"$\beta_x$")
        axs[ii].plot(tw_arcij["s"], tw_arcij["bety"], label=r"$\beta_y$")
        axs[ii].set_xlabel("s [m]")
        axs[ii].set_ylabel(r"$\beta_{x,y}$ [m]")
        axs[ii].legend()
        # axs[ii].set_ylim(-0.000001, 0.000001)
    plt.show()


def plot_xing_ip15(collider, bim=1):
    fig, axs = plt.subplots(2, figsize=(10, 15))
    for ii, i in enumerate(ips):
        tw_ipi = twiss_ip(collider, i, bim)
        axs[ii].set_title(f"IP {i}, Beam {bim}")
        axs[ii].plot(tw_ipi["s"], tw_ipi["x"], label="x", color="black", lw=3)
        axs[ii].plot(tw_ipi["s"], tw_ipi["y"], label="y", color="red", lw=3)
        axs[ii].set_xlabel("s [m]")
        axs[ii].set_ylabel(r"$co [m]$")
        axs[ii].legend()
        axs[ii].grid()
    plt.show()


def plot_xing_from_tw(tw, bim=None):
    if bim is not None:
        if bim == 1:
            tw = tw["lhcb1"]
        elif bim == 2:
            tw = tw["lhcb2"]

    fig, axs = plt.subplots(figsize=(21, 18))
    axs.grid()
    axs.set_title(f"Beam {bim}")
    axs.plot(tw["s"], tw["x"], label="x", color="black", lw=3)
    axs.plot(tw["s"], tw["y"], label="y", color="red", lw=3)
    # axs.set_ylim(-0.016, 0.016)
    axs00 = axs.twinx()
    axs00.plot(tw["s"], tw["dx"], label="Dx", color="green", lw=1, linestyle="--")
    axs00.plot(tw["s"], tw["dy"], label="Dy", color="blue", lw=1, linestyle="--")
    axs00.set_ylim(-4, 4)
    axs00.set_ylabel("$Dx,Dy \\ [m]$")
    axs01 = axs.twiny()
    axs01.set_xticks(tw[["s"], "ip.*"], tw[["name"], "ip.*"], rotation="horizontal")
    axs01.set_xlim(axs.get_xlim())
    axs.set_xlabel("$s \\ [m]$")
    axs.set_ylabel("$CO \\ [m]$")
    axs.legend(loc=2)
    axs00.legend(loc=0)
    plt.show()


def plot_xing(collider, bim=1):
    twc = collider.twiss()
    plot_xing_from_tw(twc, bim)


def plot_beta_ip15(collider, bim=1):
    fig, axs = plt.subplots(2, figsize=(20, 15))
    for ii, i in enumerate(ips):
        tw_ipi = twiss_ip(collider, i, bim)
        axs[ii].set_title(f"IP {i}, Beam {bim}")
        axs[ii].plot(
            tw_ipi["s"], tw_ipi["betx"], label=r"$\beta_x$", color="black", lw=2
        )
        axs[ii].plot(
            tw_ipi["s"], tw_ipi["bety"], label=r"$\beta_y$", color="red", lw=1.5
        )
        axs00 = axs[ii].twinx()
        axs00.plot(
            tw_ipi["s"],
            tw_ipi["dx"],
            label="Dx",
            color="green",
            lw=1,
            linestyle="--",
        )
        axs00.plot(
            tw_ipi["s"],
            tw_ipi["dy"],
            label="Dy",
            color="blue",
            lw=1,
            linestyle="--",
        )
        axs00.set_ylim(-10, 10)
        axs00.set_ylabel("Dx,Dy [m]")
        axs[ii].set_xlabel("s [m]")
        axs[ii].set_ylabel(r"$\beta_{x,y}$ [m]")
        axs[ii].legend(loc=2)
        axs[ii].grid()
        axs00.legend(loc=0)
    plt.show()


def collider_set_knobs(collider, knobs_vals, knob_vals_to_reset=None):
    if knob_vals_to_reset is not None:
        for knob_val in knob_vals_to_reset.items():
            collider.varval[knob_val[0]] = knob_val[1]

    for knob_val in knobs_vals.items():
        # print(knob_val)
        collider.varval[knob_val[0]] = knob_val[1]

    tw1 = collider.twiss()
    print("Beam1:", tw1.lhcb1[["betx", "bety", "x", "y", "px", "py"], "ip.*"])
    print("Beam2:", tw1.lhcb2[["betx", "bety", "x", "y", "px", "py"], "ip.*"])
    return tw1


# %%
config_path = "knobs_to_test.yaml"
config_knobs = ryaml.load(open(config_path, "r"))

for knob_set in config_knobs.items():
    print(knob_set)
# %%
optics_list_path = "optics_to_test.yaml"
optics_list = ryaml.load(open(optics_list_path, "r"))

for optics_set in optics_list.items():
    for optics_folder in optics_set:
        print(optics_folder)
# %%

hpre = os.getenv("HOME")
prefix = f"{hpre}/Projects/hllhc_optics/"

# optics_name = prefix + optics_list["acc-models-lhc"]["round"]["thick_150"]
optics_name = prefix + optics_list["acc-models-lhc"]["round"]["thick_580"]  # good

# optics_name = (
# prefix + "/summer_studies/collapse/opt_collapse_1000_1500_ats_500_x5hl.madx"
# )

# %%
collider = None
collider_name = "collider_hl16.json"
if os.path.exists(collider_name):
    collider = xt.Multiline.from_json(collider_name)
else:
    collider = build_collider(optics_name)
    collider.to_json(collider_name)
# %%
collider.vars.load_madx_optics_file(optics_name)
# %%
collider.lhcb1.particle_ref = xp.Particles(p0c=7000e9, q0=1, mass0=xp.PROTON_MASS_EV)
collider.lhcb2.particle_ref = xp.Particles(p0c=7000e9, q0=1, mass0=xp.PROTON_MASS_EV)
collider.build_trackers()

# %%
collider.lhcb1.twiss_default["method"] = "4d"
collider.lhcb2.twiss_default["method"] = "4d"
collider.lhcb2.twiss_default["reverse"] = True
# %%
tw0 = collider.twiss()
print("Beam1:", tw0.lhcb1[["betx", "bety"], "ip.*"])
print("Beam2:", tw0.lhcb2[["betx", "bety"], "ip.*"])

# %%
collider.lhcb1.cycle("ip3", inplace=True)
collider.lhcb2.cycle("ip3", inplace=True)
# collider_set_knobs(collider, config_knobs["xing_5_long"], config_knobs["default"])
# collider_set_knobs(collider, config_knobs["xing_5_short"], config_knobs["default"])

# %%
for knob_set in config_knobs.items():
    print(f"Testing set {knob_set[0]}")
    collider_set_knobs(collider, knob_set[1], config_knobs["default"])
    print(50 * "#")

# %%

# check arcs
for knob_set in config_knobs.items():
    print(f"Testing set {knob_set[0]} in arcs")
    collider_set_knobs(collider, knob_set[1], config_knobs["default"])
    for i, j in arcs:
        for bim in [1, 2]:
            tw_arcij = twiss_arc(collider, i, j, bim)
            assert np.isclose(np.max(tw_arcij.x), 0, atol=1e-7, rtol=0)
            assert np.isclose(np.max(tw_arcij.y), 0, atol=1e-7, rtol=0)

            assert np.isclose(np.max(tw_arcij.px), 0, atol=1e-7, rtol=0)
            assert np.isclose(np.max(tw_arcij.py), 0, atol=1e-7, rtol=0)


# check ip1
for knobs_set in [config_knobs["xing_1_short"], config_knobs["xing_1_long"]]:
    collider_set_knobs(collider, knobs_set, config_knobs["default"])
    tw_kn = collider.twiss()
    ips_tw_b1 = tw_kn.lhcb1.rows["ip.*"]
    ips_tw_b2 = tw_kn.lhcb2.rows["ip.*"]

    assert np.isclose(
        ips_tw_b1.rows["ip1"].px[0], knobs_set["on_x1"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(
        -ips_tw_b2.rows["ip1"].px[0], knobs_set["on_x1"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(ips_tw_b1.rows["ip1"].py[0], 0, atol=1e-7, rtol=0)
    assert np.isclose(ips_tw_b2.rows["ip1"].py[0], 0, atol=1e-7, rtol=0)


# %%
# check ip5
for knobs_set in [config_knobs["xing_5_short"], config_knobs["xing_5_long"]]:
    collider_set_knobs(collider, knobs_set, config_knobs["default"])
    tw_kn = collider.twiss()
    ips_tw_b1 = tw_kn.lhcb1.rows["ip.*"]
    ips_tw_b2 = tw_kn.lhcb2.rows["ip.*"]

    assert np.isclose(
        ips_tw_b1.rows["ip5"].py[0], knobs_set["on_x5"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(
        -ips_tw_b2.rows["ip5"].py[0], knobs_set["on_x5"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(ips_tw_b1.rows["ip5"].px[0], 0, atol=1e-7, rtol=0)
    assert np.isclose(ips_tw_b2.rows["ip5"].px[0], 0, atol=1e-7, rtol=0)

# %%
# check ip1 and ip5
for knobs_set in [config_knobs["xing_15_short"], config_knobs["xing_15_long"]]:
    collider_set_knobs(collider, knobs_set, config_knobs["default"])
    tw_kn = collider.twiss()
    ips_tw_b1 = tw_kn.lhcb1.rows["ip.*"]
    ips_tw_b2 = tw_kn.lhcb2.rows["ip.*"]

    assert np.isclose(
        ips_tw_b1.rows["ip5"].py[0], knobs_set["on_x5"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(
        -ips_tw_b2.rows["ip5"].py[0], knobs_set["on_x5"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(
        ips_tw_b1.rows["ip1"].px[0], knobs_set["on_x1"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(
        -ips_tw_b2.rows["ip1"].px[0], knobs_set["on_x1"] * 1e-6, atol=1e-7, rtol=0
    )
    assert np.isclose(ips_tw_b1.rows["ip1"].py[0], 0, atol=1e-7, rtol=0)
    assert np.isclose(ips_tw_b2.rows["ip1"].py[0], 0, atol=1e-7, rtol=0)
    assert np.isclose(ips_tw_b1.rows["ip5"].px[0], 0, atol=1e-7, rtol=0)
    assert np.isclose(ips_tw_b2.rows["ip5"].px[0], 0, atol=1e-7, rtol=0)

#
# %%

tw_k = collider.twiss()
# %%

plot_xing_from_tw(tw_k, bim=1)
plot_xing_from_tw(tw_k, bim=2)
# %%
plot_xing_ip15(collider, bim=1)
plot_xing_ip15(collider, bim=2)
# %%
plot_beta_ip15(collider, bim=1)
plot_beta_ip15(collider, bim=2)

# %%
plot_arcs_betas(collider, bim=1)
plot_arcs_betas(collider, bim=2)
# plot_arcs_xy(collider)

# %%
plot_arcs_xy(collider, bim=1)
plot_arcs_xy(collider, bim=2)

# %%
tw_arc12 = twiss_arc(collider, 1, 2, 1)
# %%

plot_xing(collider, bim=1)
plot_xing(collider, bim=2)

# %%
