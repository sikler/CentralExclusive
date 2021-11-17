load "../gnu/myStyleEps.gnu"

set key Left reverse samplen 1 noauto box opaque
unset colorbox

set log cb; set format cb "10^{%T}"

val_simp = `awk '{if(NR==4) print $3}' ../out/rp/polygon/res.dat`

set auto xfix
set auto yfix

set xrange [] writeback
set yrange [] writeback

set output "/dev/null"
plot "../out/rp/polygon/map.dat" w image
set output

set xrange restore
set yrange restore

set xlabel "slope"
set ylabel "intercept [pitch]"

set key title sprintf("{/Symbol c}2_{simp} = %.4f",val_simp) width -4

set xtics -1,0.05

x1 = -0.4; x2 = 0.4

set output sprintf("../eps/rp/polygon/%d.eps",j)
plot "../out/rp/polygon/map.dat"    w image, \
     "../out/rp/polygon/poly.dat" t "polygon"  w l lt 1 lc 0 lw 10, \
     "../out/rp/polygon/bands.dat" \
      u (x1):($1*x1+$2):(x2-x1):($1*(x2-x1)):($3+1) w vectors lc var nohead, \
     "../out/rp/polygon/res.dat" \
        i 0 u ($1!=0 && $2!=0 ? $1 : 1/0):2 \
                  t "centroid (poly)"          w p lt 1 pt 5, \
     "" i 1 u 1:2 t "min {/Symbol c}2 (simp)"  w p lt 4 pt 6 ps 2

