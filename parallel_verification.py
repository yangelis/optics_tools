# %%
from multiprocessing import Pool, Manager
import xtrack as xt
import xpart as xp
import xmask as xm
import xmask.lhc as xlhc

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cpymad.madx import Madx

import xtrack._temp.lhc_match as lm
import ruamel.yaml

ryaml = ruamel.yaml.YAML()
# %%


class xspoll:
    def __init__(self, collider, optics_files):
        self.collider = collider
        self.optics_files = optics_files
        self.n = len(self.optics_files)

    def start(self, func=None):
        if func is None:
            func = self.f
        manager = Manager()
        tw_xs_b1 = manager.list([pd.DataFrame for _ in range(self.n)])
        tw_xs_b2 = manager.list([pd.DataFrame for _ in range(self.n)])

        coms = [
            {
                "id": index,
                "optics_name": self.optics_files[index],
                "collider": collider,
                "txb1": tw_xs_b1,
                "txb2": tw_xs_b2,
            }
            for index in range(self.n)
        ]
        with Pool(processes=self.n) as p:
            p.map(func, coms)

        self.tw_xs = twiss_xsuite

    def f(self, mem):
        process_id = mem["id"]
        collider_to_test = mem["collider"]
        collider_to_test.vars.load_madx_optics_file(optics_name)
        tw = collider_to_test.twiss()

        mem["txb1"][process_id] = tw.lhcb1.to_pandas()
        mem["txb2"][process_id] = tw.lhcb2.to_pandas()


# %%

# config_path = "knobs_to_test.yaml"
# config_knobs = ryaml.load(open(config_path, "r"))

# for knob_set in config_knobs.items():
#     print(knob_set)
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
    collider = build_collider(optics_name)
    collider.to_json(collider_name)
# %%
optics = list(optics_list["summer_studies"]["collapse"]["round"].values())

mp = xspoll(collider, optics)
# %%
mp.start()
# %%

tw1 = xt.TwissTable(mp.tw_xs[0])
# %%
tw1[["name"], "ip.*"]

# %%
