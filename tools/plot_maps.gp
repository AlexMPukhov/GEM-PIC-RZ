# gnuplot script: 2D maps and on-axis lineout of one diagnostic step -> PDF
# usage:  gnuplot -e "step=50; dir='lwfa_witness'" tools/plot_maps.gp   (writes maps_00050.pdf)
#         gnuplot -e "step=10; dir='pwfa'; rlim=60; ezlim=0.02" tools/plot_maps.gp
#         gnuplot -e "step=1; dir='collapse_H'; rlim=15; ezlim=0.005; ions=1; nilim=20" tools/plot_maps.gp
#   rlim : half height of the maps in um (default 15)
#   ezlim: colour range of E_z in m c omega0 / e (default 0.15)
#   ions : 1 adds a map of the mobile-ion density Z n_i/n0 (default 0)
#   nilim: colour range of the ion density (default 3)
if (!exists("step")) step = 50
if (!exists("dir"))  dir  = "diag"
if (!exists("rlim")) rlim = 15
if (!exists("ezlim")) ezlim = 0.15
if (!exists("ions"))  ions  = 0
if (!exists("nilim")) nilim = 3
f   = sprintf("%s/map_%05d.txt",  dir, step)
fa  = sprintf("%s/axis_%05d.txt", dir, step)
fb  = sprintf("%s/beam_%05d.txt", dir, step)
set terminal pdfcairo enhanced size 7in,(ions ? 12 : 9) font "Helvetica,10"
set output sprintf("maps_%05d.pdf", step)
set multiplot layout (ions ? 4 : 3),1 title sprintf("GEM-PIC-RZ  %s, step %d", dir, step) noenhanced
set view map
set xlabel "{/Symbol z} [{/Symbol m}m]"; set ylabel "x [{/Symbol m}m]"
set xrange [0:*]; set yrange [-rlim:rlim]
set palette defined (-1 "blue", 0 "white", 1 "red")
set cbrange [-ezlim:ezlim]
set title "E_z  ({/Symbol f} = 0, {/Symbol p});  bunches black, trapped plasma and ionization electrons green"
splot f u 1:2:4 w pm3d notitle, \
      fb u 3:($8 >= 0 ? $1 : NaN):(0) w dots lc rgb "black" notitle, \
      fb u 3:($8 <  0 ? $1 : NaN):(0) w dots lc rgb "dark-green" notitle
set palette rgbformulae 33,13,10
set cbrange [0:4]
set title "n_e / n_0  (m = 0)"
splot f u 1:2:5 w pm3d notitle
if (ions) {
  set cbrange [0:nilim]
  set title "Z n_i / n_0  (ion blocks, m = 0)"
  splot f u 1:2:6 w pm3d notitle
}
unset view; set yrange [*:*]
set ylabel "field [m c {/Symbol w}_0 / e]"; set title "on axis"
plot fa u 1:($2/10) w l lw 1 t "E_x / 10", fa u 1:3 w l lw 1.5 t "E_z"
unset multiplot
unset output
