load ranked_0.pdb
set opaque_background
bg_color white

#color by chain
color wheat, chain A
color teal, chain B
ray 
save structure.png


as ribbon
set ribbon_width, 8
spectrum b, red_white_blue, minimum=75, maximum=90
ray
save pLDDT_structure.png



#find interface hbonds
as cartoon
color wheat, chain A
color teal, chain B

show sticks
set stick_transparency, 0.5
select interface1, chain B around 8 and (not chain B)
select interface2, chain A around 8 and (not chain A)
select interface_both, interface1 + interface2
distance test, interface1, interface2, cutoff=3.0, mode=2
hide labels
set dash_width, 5
set dash_gap, 0
orient interface_both
center interface_both
zoom interface_both
ray
save interface_hbonds.png


exit
