#!/bin/bash

for f in */*ranked_0.pdb

do

echo $f
cat > commands2plddtfigure.pml <<-EOF
load $f
orient
set_color n0, [0.051, 0.341, 0.827]
set_color n1, [0.416, 0.796, 0.945]
set_color n2, [0.996, 0.851, 0.212]
set_color n3, [0.992, 0.490, 0.302]
color n0, b < 100
color n1, b < 90
color n2, b < 70
color n3, b < 50
as cartoon
set ray_shadows, 0
set specular, off

png $f.png, width=10cm, dpi=300
EOF

pymol -Qc commands2plddtfigure.pml

done
