# Streamlit app to calculate the overlap between adjacent defocused spots
# in a STXM microscope
#
# Treats the spot as a disk or as a donut
# Disk and donut size taken from simple geometric arguments (focal length and ZP diameters
#
# Stephen Urquhart, April 19, 2026
# With a lot of CoPilot help
#

import streamlit as st
import math

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




st.title("Ptychography Overlap Calculator")
st.write("Compute the percent overlap between adjacent focus spots.")

Point_spacing = st.number_input("Step size in micron", value=0.5)
Point_spacing = Point_spacing * um

ZP_diameter = st.number_input("Zone plate diameter (micron)", value=250.0, min_value=0.0)
ZP_diameter = ZP_diameter * um

ZP_central_stop_diameter = st.number_input("Zone plate central stop diameter (micron)", value=100.0, min_value=0.0)
ZP_central_stop_diameter = ZP_central_stop_diameter * um

delta_r = st.number_input("Zone plate outer zone width (nm)", value=25.0, min_value=5.0)
delta_r = delta_r * nm

energy_eV = st.number_input("Photon Energy (eV)", value=710.0, min_value=100.0)

defocus_z_increment = st.number_input("Zone plate Z displacement for defocus (micron)", value=10.0, min_value=0.0)
Defocus_increment = defocus_z_increment * um
# -------------------------------------------------
# Zone plate parameters
# -------------------------------------------------
#ZP_diameter = 250 * um             # zone plate diameter
#delta_r = 25 * nm               # outermost zone width
#ZP_central_stop_diameter = 100 * um
#energy_eV = 710                 # photon energy

#Defocus_increment = 10 * um


#-------------------
# Math that follows
#__________________

ZP_radius = ZP_diameter / 2                # radius
CS_radius = ZP_central_stop_diameter /2

# -------------------------------------------------
# Convert photon energy to wavelength
# -------------------------------------------------
h = 6.62607015e-34               # Planck constant [J·s]
c = 2.99792458e8                # speed of light [m/s]
e = 1.602176634e-19             # elementary charge [J/eV]

wavelength = h * c / (energy_eV * e)

# -------------------------------------------------
# Zone plate focal length
# f = 2 R Δr / λ
# -------------------------------------------------
focal_length = (2 * ZP_radius * delta_r) / wavelength
defocused_focal_length = focal_length + Defocus_increment

print(f"Wavelength  = {wavelength/nm:.3f} nm")
print(f"Energy  = {energy_eV:.3f} eV")
print(f"Focal length = {focal_length/um:.0f} micron")
print(f"Defocused Focal length = {defocused_focal_length/um:.0f} micron")

# Parameters
m_outer = ZP_radius/focal_length
b_outer = -ZP_radius
m_inner = CS_radius / focal_length
b_inner = - CS_radius
defocused_outer = m_outer * defocused_focal_length + b_outer
defocused_inner = m_inner * defocused_focal_length + b_inner

#print(f"Defocused outer diameter = {defocused_outer/um:.3f} um")
#print(f"Inner donut diameter = {defocused_inner/um:.3f} um")

#disk_overlap_percent = circle_overlap_percent_largest(Point_spacing, defocused_outer, defocused_outer)

#print(f"Overlap Percentage as disk = {disk_overlap_percent:.2f}")

#donut_overlap_percent = circle_overlap_percent_largest(Point_spacing, defocused_outer, defocused_outer) - 2* circle_overlap_percent_largest(Point_spacing, defocused_outer, defocused_inner)


#print(f"Overlap Percentage as donut = {donut_overlap_percent:.2f}")


if st.button("Calculate Overlap"):
    disk_overlap_percent = circle_overlap_percent_largest(Point_spacing, defocused_outer, defocused_outer)
    donut_overlap_percent = circle_overlap_percent_largest(Point_spacing, defocused_outer,
                                                           defocused_outer) - 2 * circle_overlap_percent_largest(
        Point_spacing, defocused_outer, defocused_inner)

    st.subheader(f"Disk Overlap: {disk_overlap_percent:.2f}%")
    st.subheader(f"Donut Overlap: {donut_overlap_percent:.2f}%")
