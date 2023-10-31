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

default_tol = {
    None: 1e-8,
    "betx": 1e-6,
    "bety": 1e-6,
}  # to have no rematching w.r.t. madx


import scienceplots

# plt.rcParams["backed"] =
plt.style.use(["notebook", "science"])
plt.rcParams["legend.frameon"] = True
plt.rcParams["legend.title_fontsize"] = 20
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.labelsize"] = 20
plt.rcParams["text.usetex"] = False
plt.rcParams["xtick.labelsize"] = 20
plt.rcParams["ytick.labelsize"] = 20
plt.rcParams["legend.framealpha"] = 1

plt.rcParams["axes.titlesize"] = 22
plt.rcParams["axes.titleweight"] = "bold"

ver_lhc = 1.6
arcs = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1]]
ips = [1, 5]


def build_collider(optics_name):
    mad1 = Madx()
    mad1.call("../acc-models-lhc/lhc.seq")
    mad1.call("../acc-models-lhc/hllhc_sequence.madx")
    mad1.input("beam, sequence=lhcb1, particle=proton, energy=7000;")
    mad1.use("lhcb1")
    mad1.call(optics_name)

    mad4 = Madx()
    mad4.input("mylhcbeam=4")
    mad4.call("../acc-models-lhc/lhcb4.seq")
    mad4.call("../acc-models-lhc/hllhc_sequence.madx")
    mad4.input("beam, sequence=lhcb2, particle=proton, energy=7000;")
    mad4.use("lhcb2")
    mad4.call(optics_name)

    line1 = xt.Line.from_madx_sequence(
        mad1.sequence.lhcb1,
        allow_thick=True,
        deferred_expressions=True,
        replace_in_expr={"bv_aux": "bvaux_b1"},
    )

    line4 = xt.Line.from_madx_sequence(
        mad4.sequence.lhcb2,
        allow_thick=True,
        deferred_expressions=True,
        replace_in_expr={"bv_aux": "bvaux_b2"},
    )
    # Remove solenoids (cannot backtwiss for now)
    for ll in [line1, line4]:
        tt = ll.get_table()
        for nn in tt.rows[tt.element_type == "Solenoid"].name:
            ee_elen = ll[nn].length
            ll.element_dict[nn] = xt.Drift(length=ee_elen)

    collider = xt.Multiline(lines={"lhcb1": line1, "lhcb2": line4})
    collider.lhcb1.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7000e9)
    collider.lhcb2.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7000e9)

    collider.lhcb1.twiss_default["method"] = "4d"
    collider.lhcb2.twiss_default["method"] = "4d"
    collider.lhcb2.twiss_default["reverse"] = True

    collider.build_trackers()

    return collider


def twiss_arc(lhc: xt.Multiline, nir1: int, nir2: int, bim: int):
    tw_arc = lhc[f"lhcb{bim}"].twiss()
    return tw_arc[:, f"s.ds.r{nir1}.b{bim}":f"e.ds.l{nir2}.b{bim}"]


def twiss_ip_(lhc: xt.Multiline, nir: int, bim: int):
    tw_init = xt.TwissInit(
        element_name=f"s.ds.l{nir}.b{bim}",
        particle_ref=lhc[f"lhcb{bim}"].particle_ref,
    )
    tw_ip = lhc[f"lhcb{bim}"].twiss(
        ele_start=f"s.ds.l{nir}.b{bim}",
        ele_stop=f"e.ds.r{nir}.b{bim}",
        twiss_init=tw_init,
    )
    return tw_ip


def twiss_ip(lhc: xt.Multiline, nir: int, bim: int):
    tw_ip = lhc[f"lhcb{bim}"].twiss()
    return tw_ip[:, f"s.ds.l{nir}.b{bim}":f"e.ds.r{nir}.b{bim}"]


def build_sequence(
    mad, mylhcbeam, optics_version=1.6, ignore_cycling=False, ignore_CC=False, prefix=""
):
    # Select beam
    mad.input(f"mylhcbeam = {mylhcbeam};")

    mad.input(
        f"""
    option,-echo,-info;
    system,"mkdir temp";
    call,file="{prefix}/acc-models-lhc/lhc.seq";
    call,file="{prefix}/acc-models-lhc/toolkit/macro.madx";
    """
    )

    mad.input(
        f"""
      ! Build sequence
      option, -echo,-warn,-info;
      if (mylhcbeam==4){{
        call,file="{prefix}acc-models-lhc/lhcb4.seq";
      }} else {{
        call,file="{prefix}/acc-models-lhc/lhc.seq";
      }};
      option, -echo, warn,-info;
      """
    )

    mad.input(
        """
    l.mbh = 0.001000;
    ACSCA, HARMON := HRF400;
    """
    )

    mad.input(
        f"""
      !Install HL-LHC
      call, file=
        "{prefix}/acc-models-lhc/hllhc_sequence.madx";
      """
        """
      ! Slice nominal sequence
      exec, myslice;
      """
    )

    if not ignore_cycling:
        mad.input(
            """
        !Cycling w.r.t. to IP3 (mandatory to find closed orbit in collision in the presence of errors)
        if (mylhcbeam<3){
        seqedit, sequence=lhcb1; flatten; cycle, start=IP3; flatten; endedit;
        };
        seqedit, sequence=lhcb2; flatten; cycle, start=IP3; flatten; endedit;
        """
        )

    if not ignore_CC:
        mad.input(
            f"""
        ! Install crab cavities (they are off)
        call, file='{prefix}/acc-models-lhc/toolkit/enable_crabcavities.madx';
        on_crab1 = 0;
        on_crab5 = 0;
        """
        )

    mad.input(
        """
        ! Set twiss formats for MAD-X parts (macro from opt. toolkit)
        exec, twiss_opt;
        """
    )


def apply_optics(mad, optics_file):
    mad.call(optics_file)
    # A knob redefinition
    mad.input("on_alice := on_alice_normalized * 7000./nrj;")
    mad.input("on_lhcb := on_lhcb_normalized * 7000./nrj;")


def build_collider_from_mad(optics_filename, ver_lhc=1.6, prefix=""):
    # Make mad environment
    links = {"acc-models-lhc": "../"}
    xm.make_mad_environment(links=links)

    # Start mad
    mad_b1b2 = Madx(command_log="mad_collider.log")
    mad_b4 = Madx(command_log="mad_b4.log")

    # Build sequences
    build_sequence(
        mad_b1b2,
        mylhcbeam=1,
        optics_version=ver_lhc,
        prefix=prefix,
        ignore_cycling=True,
    )
    build_sequence(
        mad_b4, mylhcbeam=4, optics_version=ver_lhc, prefix=prefix, ignore_cycling=True
    )

    # Apply optics (only for b1b2, b4 will be generated from b1b2)
    apply_optics(mad_b1b2, optics_file=optics_filename)

    beam_config = {
        "lhcb1": {"beam_energy_tot": 7000},
        "lhcb2": {"beam_energy_tot": 7000},
    }
    # Build xsuite collider
    collider = xlhc.build_xsuite_collider(
        sequence_b1=mad_b1b2.sequence.lhcb1,
        sequence_b2=mad_b1b2.sequence.lhcb2,
        sequence_b4=mad_b4.sequence.lhcb2,
        beam_config=beam_config,
        enable_imperfections=False,
        enable_knob_synthesis=True,
        rename_coupling_knobs=True,
        # pars_for_imperfections=None,
        # ver_lhc_run=None,
        ver_hllhc_optics=ver_lhc,
    )

    # Return collider
    return collider


def detuning_terms(collider, mo_list, bim=1):
    results = {
        "dqxdex": [],
        "dqxdey": [],
        "dqydey": [],
    }

    for mo in mo_list:
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

        results["dqxdex"].append(dqxdex)
        results["dqxdey"].append(dqxdey)
        results["dqydey"].append(dqydey)
    collider[f"lhcb{bim}"].varval[f"i_oct_b{bim}"] = 0
    return results
