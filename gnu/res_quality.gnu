load "myStyleEps.gnu"

set datafile missing "0"

f(x) = a*x
a = 1

data = "<res_quality.awk \
         ../out/res/rp_phi_pi-low_TB.his \
         ../out/res/rp_phi_pi-low_BT.his \
         ../out/res/rp_phi_pi-low_TT.his \
         ../out/res/rp_phi_pi-low_BB.his" 

set label 9 at graph 0.05,0.8 left

set style data p

ps = 1

#
set size ratio -1
set xrange [0:12]
set yrange [0:12]

# TB vs BT
set output "../eps/res/rp_phi_qual_a.eps"
set xlabel "{/Symbol s}(TB)"
set ylabel "{/Symbol s}(BT)"

fit f(x) data u 1:2:5:6 xyerrors via a
set label 9 sprintf("a = %.3f +- %.3f",a,a_err)
plot       "" u 1:2 pt 6 ps ps, f(x) lt 3

# BB vs TT
set xrange [0:18]
set yrange [0:18]

set output "../eps/res/rp_phi_qual_b.eps"
set xlabel "{/Symbol s}(TT)"
set ylabel "{/Symbol s}(BB)"

fit f(x) "" u 3:4:7:8 xyerrors via a
set label 9 sprintf("a = %.3f +- %.3f",a,a_err)
plot     "" u 3:4 pt 6 ps ps, f(x) lt 3

#
set xrange [0:6]
set yrange [0:6]

# TB vs TT
set output "../eps/res/rp_phi_qual_c.eps"
set xlabel "{/Symbol s}(TT)"
set ylabel "{/Symbol s}(TB)"

fit f(x) "" u 3:1:7:5 xyerrors via a
set label 9 sprintf("a = %.3f +- %.3f",a,a_err)
plot     "" u 3:1 pt 6 ps ps, f(x) lt 3

# BT vs TT
set output "../eps/res/rp_phi_qual_d.eps"
set xlabel "{/Symbol s}(TT)"
set ylabel "{/Symbol s}(BT)"

fit f(x) "" u 3:2:7:6 xyerrors via a
set label 9 sprintf("a = %.3f +- %.3f",a,a_err)
plot     "" u 3:2 pt 6 ps ps, f(x) lt 3
