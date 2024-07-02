
#This is an example Script for proceessing the docked protein-ligand complex. Download Interaction folder complex structure, and run this
#in pymol. You still need to manually adding features like distance of H bond, or add H atoms to certain residues with h_add
# commond etc. It is not a magic script to handle everything but good enough to finish 70% of the work to produce a publication
# quality image. You can turn on/off certain lines by adding or removing #.

load output.pdb, tmp
#load Docked1.pdb, docked
#select the surrounding then set it to an object 
# find 5 angstroms around ligand
select resn UNL
select res, byres sele around 5
#set the ligand to an lig object
select lig, resn UNL
# Select all to hide
select all
# hide everythong
hide
# show the ligand and surroundings
show sticks, lig|res
# hide the H atom on the sticks
hide sticks, h.
# zoon on lig and res, background faded
zoom lig|res
# Color ligand carbon atom yellow
color yellow, lig&name C*
# Color ligand carbon atom yellow
color cyan, res&name C*
# Set label and font size
# note you still have to click "residues" form L dropdown menue to show the label
set label_size, 28
set label_font_id, 8
set label_color, grey

##The next secotion is to display HYDOPHOBIC interacton (but it will not finish all the work, mainly make the psudoatom for your automatically).
# First generate pesudo atoms/ please adjust accordingly/ Here we show the pseudoatom for Gatekeeper of this kinase
select /tmp//A/PHE`241/CD2 | /tmp//A/PHE`241/CE2 | /tmp//A/PHE`241/CZ | /tmp//A/PHE`241/CE1 | /tmp//A/PHE`241/CD1 | /tmp//A/PHE`241/CG ##Change to your own site numbers!
pseudoatom pi_center, sele
#select res, pi_center
show sphere, pi_center
set sphere_scale, 0.25
set sphere_color, cyan

# We also define some hydrophobic interaction style
set sphere_scale, 0.35
set sphere_transparency, 0.3
set sphere_color, cyan
# Define some hydrogen bond. interaciton style
set dash_length, 0.3
set dash_radius, 0.08
set dash_gap, 0.3
# set backgroud white
bg_color white
show cartoon
set cartoon_transparency, 0.3
bg_color white
set ray_trace_mode, 1
set_view (\
     0.741649985,    0.022153847,   -0.670419037,\
    -0.523934960,    0.643225133,   -0.558340371,\
     0.418860048,    0.765345991,    0.488656282,\
     0.000000000,    0.000000000,  -60.401882172,\
    11.944400787,   53.451717377,   28.071718216,\
    54.113571167,   66.690200806,  -20.000000000 )
png binding_pocket.png
save pre_pymol.pse
