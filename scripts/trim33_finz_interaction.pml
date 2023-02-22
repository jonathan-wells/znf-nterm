reinit

/* util.performance(0) */
bg white
set specular, 0.1
set ray_trace_mode, 1
set antialias, 3

load /Users/jonwells/Projects/feschottelab/znf-nterm/data/structures/TRIM33_FINZ.pdb
split_chains
set_name TRIM33_FINZ_A, TRIM33_A
set_name TRIM33_FINZ_B, TRIM33_B
set_name TRIM33_FINZ_C, FINZ
delete TRIM33_FINZ

/* set base colours */
color gray70, TRIM33_A
color gray70, TRIM33_B
set surface_color, gray60, TRIM33_A
set surface_color, gray60, TRIM33_B
color smudge, FINZ
set surface_color, smudge, FINZ
set transparency, 0.70
show surface

/* show contacted residues */
contacts FINZ, TRIM33_A, result='A_contacts', bigcutoff=4.0
contacts FINZ, TRIM33_B, result='B_contacts', bigcutoff=4.0
select finz_res, FINZ and res 28+29+32
select a_res, TRIM33_A and res 53

/* set stick_color, white, (finz_res or a_res) and elem C */
/* set sphere_color, white, (finz_res or a_res) and elem C */
show sticks, (finz_res or a_res) and (sidechain or name CA)
show spheres, (finz_res or a_res) and (sidechain or name CA)

set stick_radius, .07
set sphere_scale, .22
set sphere_scale, .13, elem H
set stick_quality, 50
set sphere_quality, 4
/* color gray85, elem C */
color red, elem O
color slate, elem N
color gray98, elem H


