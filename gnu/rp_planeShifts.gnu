load "myStyleEps.gnu"

#
set term push

set term post eps enh color dashed dl 2 "Helvetica" 30 size 3*5in,3*4.5in

set output "../eps/rp/planeShifts_line.eps"

set multiplot layout 4,4 margins 0.08,0.92,0.06,0.99 spacing 0.06

set xlabel "Relative shift [pitch]"
#set ylabel "Joint goodness-of-fit [10^6]" offset 2,0

set ytics  nomirror
set y2tics nomirror

set style data l

do for [i=0:15] {
 set label 9 lab(i) at graph 0.95,0.9 right
 plot [-1:1] \
  "../out/rp/2part/planeShifts_chi2_line.dat" \
     i i u ($3==1 ? $1 : 1/0):($2*1e-6) ls 1 dt 1 axes x1y1, \
  "" i i u ($3==2 ? $1 : 1/0):($2*1e-6) ls 2 dt 1 axes x1y1, \
  "" i i u ($3==3 ? $1 : 1/0):($2*1e-6) ls 3 dt 1 axes x1y1, \
  "../out/rp/2part/planeShifts_zero_line.dat" \
     i i u ($3==1 ? $1 : 1/0):($2*1e-7) ls 1 dt 3 axes x1y2, \
  "" i i u ($3==2 ? $1 : 1/0):($2*1e-7) ls 2 dt 3 axes x1y2, \
  "" i i u ($3==3 ? $1 : 1/0):($2*1e-7) ls 3 dt 3 axes x1y2, \
}

unset multiplot
set output

unset y2tics
set ytics mirror

set style data p

set xlabel
unset label 9

set term pop

#
set term post eps enh color dashed dl 2 "Helvetica" 30 size 2*5in,4.5in

set key bottom left box opaque samplen 0 width -1 \
    maxrows 3 # title "joint {/Symbol c}^2 | count 0s"

do for [i=0:15] {
 set xtics add (lab(i) i) offset 0,-2 rotate by 90
}

set ylabel "Relative shift [pitch]" offset 0,0

do for [i=0:16:2] {
  set arrow from i-0.5,graph 0 rto 0,graph 1 lt 1 dt 2 lc 0 nohead
}

dx = 0.05

set output "../eps/rp/planeShifts.eps"
plot [0-0.5:16-0.5] \
 0 lt 0, \
 "../out/rp/2part/planeShifts_chi2.dat" \
    u ($0-dx):2 t "plane 2" ls 1 ps 2, \
 "" u ($0   ):3 t "plane 3" ls 2 ps 2, \
 "" u ($0+dx):4 t "plane 4" ls 3 ps 2, \
 "../out/rp/0part/planeShifts_chi2.dat" \
    u ($0-dx):2 t "plane 2" ls 1 ps 1, \
 "" u ($0   ):3 t "plane 3" ls 2 ps 1, \
 "" u ($0+dx):4 t "plane 4" ls 3 ps 1
# \
#"../out/rp/planeShifts_9_chi2.dat" \
#   u ($0-dx):2 t "plane 2" ls 1 ps 4, \
#"" u ($0   ):3 t "plane 3" ls 2 ps 4, \
#"" u ($0+dx):4 t "plane 4" ls 3 ps 4
