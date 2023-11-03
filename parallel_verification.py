# %%
from misc import *
from multiprocessing import Pool, Manager
from pprint import pprint
import ruamel.yaml

ryaml = ruamel.yaml.YAML()


# %%
def collider_set_knobs(collider, knobs_vals, knob_vals_to_reset=None):
    if knob_vals_to_reset is not None:
        for knob_val in knob_vals_to_reset.items():
            collider.varval[knob_val[0]] = knob_val[1]

    for knob_val in knobs_vals.items():
        # print(knob_val)
        collider.varval[knob_val[0]] = knob_val[1]

    tw1 = collider.twiss()
    # print("Beam1:", tw1.lhcb1[["betx", "bety", "x", "y", "px", "py"], "ip.*"])
    # print("Beam2:", tw1.lhcb2[["betx", "bety", "x", "y", "px", "py"], "ip.*"])
    return tw1


class xspoll:
    def __init__(self, collider, optics_files, knob_configs):
        self.collider = collider
        self.optics_files = optics_files
        self.n = len(self.optics_files)
        self.knob_configs = knob_configs

    def start(self, func=None):
        if func is None:
            func = self.f
        manager = Manager()
        tw_xs_b1 = manager.list([pd.DataFrame for _ in range(self.n)])
        tw_xs_b2 = manager.list([pd.DataFrame for _ in range(self.n)])

        coms = [
            {
                "id": index,
                "config": self.knob_configs,
                "optics_name": self.optics_files[index],
                "collider": collider,
                "txb1": tw_xs_b1,
                "txb2": tw_xs_b2,
            }
            for index in range(self.n)
        ]
        with Pool(processes=self.n) as p:
            p.map(func, coms)

        self.twb1 = tw_xs_b1
        self.twb2 = tw_xs_b2

    def f(self, mem):
        process_id = mem["id"]
        collider_to_test = mem["collider"]
        optics_name = mem["optics_name"]
        print(f"Checking optics file {optics_name}")

        collider_to_test.vars.load_madx_optics_file(optics_name)
        tw = collider_to_test.twiss()

        for knobs_set in [config_knobs["xing_1_short"], config_knobs["xing_1_long"]]:
            collider_set_knobs(collider_to_test, knobs_set, config_knobs["default"])
            tw_kn = collider_to_test.twiss()
            ips_tw_b1 = tw_kn.lhcb1.rows["ip.*"]
            ips_tw_b2 = tw_kn.lhcb2.rows["ip.*"]

            if not np.isclose(
                ips_tw_b1.rows["ip1"].px[0],
                knobs_set["on_x1"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b1 px= {ips_tw_b1.rows["ip1"].px[0]}, {optics_name=}')

            if not np.isclose(
                -ips_tw_b2.rows["ip1"].px[0],
                knobs_set["on_x1"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b2 px= {ips_tw_b2.rows["ip1"].px[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b1.rows["ip1"].py[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b1 py= {ips_tw_b1.rows["ip1"].py[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b2.rows["ip1"].py[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b2 py= {ips_tw_b2.rows["ip1"].py[0]}, {optics_name=}')

        for knobs_set in [config_knobs["xing_5_short"], config_knobs["xing_5_long"]]:
            collider_set_knobs(collider_to_test, knobs_set, config_knobs["default"])
            tw_kn = collider_to_test.twiss()
            ips_tw_b1 = tw_kn.lhcb1.rows["ip.*"]
            ips_tw_b2 = tw_kn.lhcb2.rows["ip.*"]

            if not np.isclose(
                ips_tw_b1.rows["ip5"].py[0],
                knobs_set["on_x5"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b1 py= {ips_tw_b1.rows["ip5"].py[0]}, {optics_name=}')

            if not np.isclose(
                -ips_tw_b2.rows["ip5"].py[0],
                knobs_set["on_x5"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b2 py= {ips_tw_b2.rows["ip5"].py[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b1.rows["ip5"].px[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b1 px= {ips_tw_b1.rows["ip5"].px[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b2.rows["ip5"].px[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b2 px= {ips_tw_b2.rows["ip5"].px[0]}, {optics_name=}')

        for knobs_set in [config_knobs["xing_15_short"], config_knobs["xing_15_long"]]:
            collider_set_knobs(collider_to_test, knobs_set, config_knobs["default"])
            tw_kn = collider_to_test.twiss()
            ips_tw_b1 = tw_kn.lhcb1.rows["ip.*"]
            ips_tw_b2 = tw_kn.lhcb2.rows["ip.*"]

            if not np.isclose(
                ips_tw_b1.rows["ip5"].py[0],
                knobs_set["on_x5"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b1 py= {ips_tw_b1.rows["ip5"].py[0]}, {optics_name=}')

            if not np.isclose(
                -ips_tw_b2.rows["ip5"].py[0],
                knobs_set["on_x5"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b2 py= {ips_tw_b2.rows["ip5"].py[0]}, {optics_name=}')

            if not np.isclose(
                ips_tw_b1.rows["ip1"].px[0],
                knobs_set["on_x1"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b1 px= {ips_tw_b1.rows["ip1"].px[0]}, {optics_name=}')

            if not np.isclose(
                -ips_tw_b2.rows["ip1"].px[0],
                knobs_set["on_x1"] * 1e-6,
                atol=1e-7,
                rtol=0,
            ):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b2 px= {ips_tw_b2.rows["ip1"].px[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b1.rows["ip1"].py[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b1 py= {ips_tw_b1.rows["ip1"].py[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b2.rows["ip1"].py[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip1 b2 py= {ips_tw_b2.rows["ip1"].py[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b1.rows["ip5"].px[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b1 py= {ips_tw_b1.rows["ip5"].px[0]}, {optics_name=}')

            if not np.isclose(ips_tw_b2.rows["ip5"].px[0], 0, atol=1e-7, rtol=0):
                print(f"Failed with config {knobs_set}")
                print(f'ip5 b2 py= {ips_tw_b2.rows["ip5"].px[0]}, {optics_name=}')

        # check all xing

        knobs_set = config_knobs["xing_all"]
        collider_set_knobs(collider_to_test, knobs_set, config_knobs["default"])
        tw_kn = collider_to_test.twiss()
        ips_tw_b1 = tw_kn.lhcb1.rows["ip.*"]
        ips_tw_b2 = tw_kn.lhcb2.rows["ip.*"]

        if not np.isclose(
            ips_tw_b1.rows["ip5"].py[0],
            knobs_set["on_x5"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip5 b1 py= {ips_tw_b1.rows["ip5"].py[0]}, {optics_name=}')

        if not np.isclose(
            -ips_tw_b2.rows["ip5"].py[0],
            knobs_set["on_x5"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip5 b2 py= {ips_tw_b2.rows["ip5"].py[0]}, {optics_name=}')

        if not np.isclose(
            ips_tw_b1.rows["ip1"].px[0],
            knobs_set["on_x1"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip1 b1 px= {ips_tw_b1.rows["ip1"].px[0]}, {optics_name=}')

        if not np.isclose(
            -ips_tw_b2.rows["ip1"].px[0],
            knobs_set["on_x1"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip1 b2 px= {ips_tw_b2.rows["ip1"].px[0]}, {optics_name=}')

        if not np.isclose(
            ips_tw_b1.rows["ip2"].py[0],
            knobs_set["on_x2"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip2 b1 px= {ips_tw_b1.rows["ip2"].py[0]}, {optics_name=}')

        if not np.isclose(
            -ips_tw_b2.rows["ip2"].py[0],
            knobs_set["on_x2"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip2 b2 px= {ips_tw_b2.rows["ip2"].py[0]}, {optics_name=}')

        if not np.isclose(
            ips_tw_b1.rows["ip8"].py[0],
            knobs_set["on_x8v"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip8 b1 px= {ips_tw_b1.rows["ip8"].py[0]}, {optics_name=}')

        if not np.isclose(
            -ips_tw_b2.rows["ip8"].py[0],
            knobs_set["on_x8v"] * 1e-6,
            atol=1e-7,
            rtol=0,
        ):
            print(f"Failed with config {knobs_set}")
            print(f'ip8 b2 px= {ips_tw_b2.rows["ip8"].py[0]}, {optics_name=}')

        mem["txb1"][process_id] = tw.lhcb1.to_pandas()
        mem["txb2"][process_id] = tw.lhcb2.to_pandas()


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
collider = None
collider_name = "collider_hl16.json"
if os.path.exists(collider_name):
    collider = xt.Multiline.from_json(collider_name)
else:
    collider = build_collider(optics_list["acc-models-lhc"]["round"]["thick_150_1500"])
    collider.to_json(collider_name)
# %%
optics = list(optics_list["summer_studies"]["collapse"]["round"].values())
pprint(optics)
mp = xspoll(collider, optics, config_knobs)
# %%
mp.start()

# %%
optics = list(optics_list["summer_studies"]["collapse"]["flat"].values())
pprint(optics)
mp = xspoll(collider, optics, config_knobs)
# %%
mp.start()
# %%
optics = list(optics_list["acc-models-lhc"]["round"].values())
pprint(optics)
mp = xspoll(collider, optics, config_knobs)
# %%
mp.start()
# %%
optics = list(optics_list["acc-models-lhc"]["flat"].values())
pprint(optics)
mp = xspoll(collider, optics, config_knobs)
# %%
mp.start()
# %%
optics = list(optics_list["acc-models-lhc"]["levelling"].values())
pprint(optics)
mp = xspoll(collider, optics, config_knobs)
# %%
mp.start()
# %%
