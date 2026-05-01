# Streamlit app to calculate the overlap between adjacent defocused spots
# in a STXM microscope
#
# Treats the spot as a disk or as a donut
# Disk and donut size taken from simple geometric arguments (focal length and ZP diameters)
#
# Stephen Urquhart
# With a lot of help from CoPilot (initially) and then Claude AI (May 1, 2026)
# Started April 19, 2026
# Updated May 1, 2026
#
# MIT License
#
# Copyright (c) 2026 Stephen Urquhart
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import streamlit as st
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

st.write("Python:", sys.executable)

um = 1e-6
nm = 1e-9

def circle_overlap_percent_largest(d, r1, r2):
    """
    Calculate percent overlap relative to the largest circle.

    Parameters
    ----------
    d : float
        Distance between circle centers
    r1 : float
        Radius of circle 1
    r2 : float
        Radius of circle 2

    Returns
    -------
    float
        Percent of the largest circle that is overlapped
    """

    # No overlap
    if d >= (r1 + r2):
        return 0.0

    # One circle completely inside the other
    if d <= abs(r1 - r2):
        return 100.0 if min(r1, r2) == max(r1, r2) else \
               100.0 * (math.pi * min(r1, r2)**2) / (math.pi * max(r1, r2)**2)

    # Partial overlap
    overlap = (
        r1**2 * math.acos((d**2 + r1**2 - r2**2) / (2 * d * r1)) +
        r2**2 * math.acos((d**2 + r2**2 - r1**2) / (2 * d * r2)) -
        0.5 * math.sqrt(
            (-d + r1 + r2) *
            ( d + r1 - r2) *
            ( d - r1 + r2) *
            ( d + r1 + r2)
        )
    )

    largest_area = math.pi * max(r1, r2)**2
    return 100.0 * overlap / largest_area


_C1  = np.array([0.25, 0.50, 0.85])   # blue  — spot 1
_C2  = np.array([0.18, 0.70, 0.45])   # green — spot 2
_COV = np.array([0.93, 0.35, 0.08])   # orange — overlap
_OV_ALPHA = 0.30                        # how opaque the overlap tint is


def plot_overlap(r, d, r_inner=0):
    """Rasterise two equal annuli (or disks) separated by d.

    Spot 1 is blue, spot 2 is green; the overlap is a semi-transparent orange
    blended over the average of the two spot colours so circle edges stay
    visible through the overlap region.
    """
    margin = max(0.2 * r, 0.5 * d)
    xs = np.linspace(-r - d / 2 - margin, r + d / 2 + margin, 700)
    ys = np.linspace(-r - margin, r + margin, 700)
    X, Y = np.meshgrid(xs, ys)

    in_c1 = (X + d / 2) ** 2 + Y ** 2 <= r ** 2
    in_c2 = (X - d / 2) ** 2 + Y ** 2 <= r ** 2
    if r_inner > 0:
        in_c1 &= (X + d / 2) ** 2 + Y ** 2 >= r_inner ** 2
        in_c2 &= (X - d / 2) ** 2 + Y ** 2 >= r_inner ** 2

    img = np.ones((*X.shape, 3))            # RGB, white background
    img[in_c1 & ~in_c2] = _C1
    img[in_c2 & ~in_c1] = _C2
    # Overlap: blend orange over the average of the two spot colours
    img[in_c1 & in_c2] = _OV_ALPHA * _COV + (1 - _OV_ALPHA) * (_C1 + _C2) / 2

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.imshow(img,
              extent=[xs[0] / um, xs[-1] / um, ys[0] / um, ys[-1] / um],
              origin='lower', aspect='equal')

    for cx in [-d / 2 / um, d / 2 / um]:
        ax.add_patch(mpatches.Circle((cx, 0), r / um,
                                     fill=False, edgecolor='#111111', linewidth=1.5, zorder=3))
        if r_inner > 0:
            ax.add_patch(mpatches.Circle((cx, 0), r_inner / um,
                                         fill=False, edgecolor='#111111', linewidth=1.5, zorder=3))

    ov_display = _OV_ALPHA * _COV + (1 - _OV_ALPHA) * (_C1 + _C2) / 2
    ax.legend(handles=[
        mpatches.Patch(color=_C1,        label='Spot 1'),
        mpatches.Patch(color=_C2,        label='Spot 2'),
        mpatches.Patch(color=ov_display, label='Overlap'),
    ], loc='upper right', fontsize=8)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title(f"r = {r/um:.2f} µm,  step = {d/um:.2f} µm")
    return fig


st.title("Ptychography Overlap Calculator")
st.write("Compute the percent overlap between adjacent focus spots.")

with st.sidebar:
    st.header("Scan Parameters")
    Point_spacing = st.number_input("Step size (micron)", value=0.5, min_value=0.0) * um
    energy_eV = st.number_input("Photon Energy (eV)", value=710.0, min_value=100.0)

    st.header("Zone Plate Parameters")
    ZP_diameter = st.number_input("Diameter (micron)", value=250.0, min_value=0.0) * um
    ZP_central_stop_diameter = st.number_input("Central stop diameter (micron)", value=100.0, min_value=0.0) * um
    delta_r = st.number_input("Outer zone width (nm)", value=25.0, min_value=5.0) * nm



# Shared physics calculations
ZP_radius = ZP_diameter / 2
CS_radius = ZP_central_stop_diameter / 2

h = 6.62607015e-34
c = 2.99792458e8
e = 1.602176634e-19
wavelength = h * c / (energy_eV * e)
focal_length = (2 * ZP_radius * delta_r) / wavelength

st.info(f"Wavelength: {wavelength/nm:.3f} nm  |  Focal length: {focal_length/um:.0f} µm")

col1, col2 = st.columns(2)

with col1:
    st.subheader("Method 1: Z-displacement")
    defocus_z = st.number_input("Z displacement (micron)", value=10.0, min_value=0.0)
    Defocus_increment = defocus_z * um
    if st.button("Calculate", key="calc_z"):
        defocused_focal_length = focal_length + Defocus_increment
        defocused_outer = (ZP_radius / focal_length) * defocused_focal_length - ZP_radius
        defocused_inner = (CS_radius / focal_length) * defocused_focal_length - CS_radius
        st.session_state.results = {
            "method": 1,
            "r_outer": defocused_outer,
            "r_inner": defocused_inner,
            "disk_overlap": circle_overlap_percent_largest(Point_spacing, defocused_outer, defocused_outer),
            "donut_overlap": (
                circle_overlap_percent_largest(Point_spacing, defocused_outer, defocused_outer)
                - 2 * circle_overlap_percent_largest(Point_spacing, defocused_outer, defocused_inner)
            ),
            "spacing": Point_spacing,
        }

with col2:
    st.subheader("Method 2: Direct Spot Size")
    defocus_spot_size = st.number_input("Defocus spot diameter (micron)", value=2.5, min_value=0.0) * um
    if st.button("Calculate", key="calc_spot"):
        spot_inner = defocus_spot_size * CS_radius / ZP_radius
        st.session_state.results = {
            "method": 2,
            "r_outer": defocus_spot_size,
            "r_inner": spot_inner,
            "z_displacement": defocus_spot_size * focal_length / ZP_radius,
            "disk_overlap": circle_overlap_percent_largest(Point_spacing, defocus_spot_size, defocus_spot_size),
            "donut_overlap": (
                circle_overlap_percent_largest(Point_spacing, defocus_spot_size, defocus_spot_size)
                - 2 * circle_overlap_percent_largest(Point_spacing, defocus_spot_size, spot_inner)
            ),
            "spacing": Point_spacing,
        }

if "results" in st.session_state:
    res = st.session_state.results
    st.divider()

    if res["method"] == 1:
        st.subheader("Results — Method 1: Z-displacement")
        m1, m2, m3 = st.columns(3)
        m1.metric("Defocused diameter", f"{res['r_outer']/um:.2f} µm")
        m2.metric("Disk overlap", f"{res['disk_overlap']:.2f}%")
        m3.metric("Donut overlap", f"{res['donut_overlap']:.2f}%")
    else:
        st.subheader("Results — Method 2: Direct Spot Size")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Defocused diameter", f"{res['r_outer']/um:.2f} µm")
        m2.metric("Required Z displacement", f"{res['z_displacement']/um:.2f} µm")
        m3.metric("Disk overlap", f"{res['disk_overlap']:.2f}%")
        m4.metric("Donut overlap", f"{res['donut_overlap']:.2f}%")

    pc1, pc2 = st.columns(2)
    with pc1:
        st.caption("Disk")
        st.pyplot(plot_overlap(res["r_outer"], res["spacing"]))
    with pc2:
        st.caption("Donut")
        st.pyplot(plot_overlap(res["r_outer"], res["spacing"], r_inner=res["r_inner"]))


