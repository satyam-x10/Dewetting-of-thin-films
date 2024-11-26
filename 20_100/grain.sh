rm -rf output45 
mkdir output45 

for i in `seq 0 2000 200000`
do

gnuplot << EOF
set terminal pngcairo size 800,800 enhanced font "Arial,22" linewidth 5
set key samplen 0 spacing 0.9 
set border lw 0.6
set xrange [0:512]
unset key
unset xtic
unset ytic
set yrange[0:512]

unset colorbox
set output 'output45/$i.png'


set pm3d map 
splot "output/grain/2d_grain$i.dat" u 1:2:5
EOF

done
