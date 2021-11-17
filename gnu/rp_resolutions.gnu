set fit quiet log "/dev/null"

#
set samples 1000

f(x) = a*exp(-((x-m)/s)**2/2)
w(x) = (x > 1 ? sqrt(x) : 1)

do for [topo in "TB BT"] {
 print " fitting ",topo

 a = 10; m = 0.01; s = 0.1

 set print sprintf("../out/rp/0part/sum2px_vs_pred_%s.par",topo)

do for [ix=95:650:10] {
  x = ix*1e-4

  fit f(x) \
   sprintf("<awk '{if($1==%f) print}' \
                    ../out/rp/0part/sum2px_vs_pred_%s.his",x,topo) \
   u 2:3:(w($4)) zerror via a,m,s

  set label 8 sprintf("%.3f",x) at graph 0.1,0.9
# plot "" u 2:3 w his, f(x) lt 3 #; pause -1

  print x, m,m_err, abs(s),s_err, FIT_WSSR, FIT_NDF
 }

 set print
}
unset label 8

#
load "myStyleEps.gnu"
set size square

c = 1e-3

set term post eps enh color dashed dl 2 "Helvetica" 30 size 2*5in,2*4.5in
set output "../eps/rp/rp_resolution_y.eps"
set multiplot layout 2,2 margins 0.08,0.98, 0.15,0.98

f(x) = g(x)*h(x)

g(x) = a*exp(-((x-m)/s)**2/2)
h(x) = 1 - b*exp(-((x-n)/z)**2/2)

w(x,e) = (x < -0.02 || x > 0.01 ? e : 1e+30)

set xlabel "{/Symbol S}p_y [GeV]"
set ylabel "counts [10^3]"

set label 5 at graph 0.05,0.9 left

#
do for [topo in "TB BT"] {
 print " pars ",topo
 set label 5 topo

 a = 1e+4; m = 0.01; s = 0.02
 b = 0.1 ; n = 0.01; z = 0.02

 fit [-0.06:0.06] f(x) \
   sprintf("../out/rp/0part/sum2py_%s.his",topo) \
   u 1:(c*$2):(w($1,$3)) zerror via a,m,s

 fit [-0.06:0.06] f(x) \
   sprintf("../out/rp/0part/sum2py_%s.his",topo) \
   u 1:(c*$2):3          zerror via a,m,s, b,n,z

 set label 6 \
  sprintf("m_y = %.1f MeV\n{/Symbol s}_y = %.1f MeV",m*1e+3,s*1e+3) \
  at graph 0.95,0.5 right

 plot [-0.15:0.15] "" u 1:(c*$2) w his, g(x) lt 3
}

#
do for [topo in "TT BB"] {
 print " pars ",topo
 set label 5 topo

 a = 1.0; m = 0.01; s = 0.02
 b = 0.1; d = 1
 f(x) = a*exp(-((x-m)/s)**2/2) + b*x + d

 fit [-0.15:0.15] f(x) \
   sprintf("../out/rp/2part/sum2py_%s.his",topo) \
   u 1:(c*$2):(w($1,$3)) zerror via a,d

 fit [-0.15:0.15] f(x) \
   sprintf("../out/rp/2part/sum2py_%s.his",topo) \
   u 1:(c*$2):3          zerror via a,m,s, b,d

 set label 6 \
  sprintf("m_y = %.1f MeV\n{/Symbol s}_y = %.1f MeV",m*1e+3,s*1e+3) \
  at graph 0.95,0.5 right

 plot [-0.15:0.15] sprintf("../out/rp/2part/sum2py_%s.his",topo) \
     u 1:(c*$2) w his, f(x) lt 3
}


unset multiplot

unset label 5

set size nosquare

set key top left

#
set term post eps enh color dashed dl 2 "Helvetica" 30 size 5in,4.5in
set output "../eps/rp/rp_resolution_x.eps"

set xtics 0,0.02
set ytics 0,0.02

set xlabel "{/Symbol s}_{{/Symbol S}p_x} from position resolution [GeV]"
set ylabel "Measured {/Symbol s}_{{/Symbol S}p_x} [GeV]" offset 2,0

f(x) = sqrt(s**2 + x**2)
s = 40e-3

# FIXME
fit [0.015:0.060] f(x) \
    "<cat ../out/rp/0part/sum2px_vs_pred_*.par" u 1:4:5 zerror via s

set label 6 sprintf("{/Symbol s}_0 = %.1f MeV",s*1e+3) \
    at graph 0.95,0.1 right

#set ytics 0.04,0.005
#set xtics 0,0.01

#set ylabel offset 2.75,0

set ytics  nomirror
set y2tics nomirror

set rmargin 6

set y2range [0:]

set y2label "counts" offset -1,0

c = 1e-3
plot [0:0.08][0:0.08] \
  "../out/rp/0part/sum2px_vs_pred_TB.par" u 1:4:5 t "TB" w e ls 1 ps 1.5, \
  "../out/rp/0part/sum2px_vs_pred_BT.par" u 1:4:5 t "BT" w e ls 3 ps 1.5, \
  f(x) ls 2 dt 1, s lt 0, \
  "../out/rp/0part/sum2px_pred_TB.his" u 1:(c*$2) axes x1y2 t "TB" w l ls 1, \
  "../out/rp/0part/sum2px_pred_BT.his" u 1:(c*$2) axes x1y2 t "BT" w l ls 3, \
  x ls 2

#
unset label 6
set ytics -0.02,0.01 mirror
unset y2tics
unset y2label
set auto y

set ylabel "Mean {/Symbol S}p_x [GeV]"

set output "../eps/rp/rp_mean_x.eps"
plot [][-0.020:0.020] \
 "../out/rp/0part/sum2px_vs_pred_TB.par" u 1:2:3 t "TB" w e ls 1, \
 "../out/rp/0part/sum2px_vs_pred_BT.par" u 1:2:3 t "BT" w e ls 3, \
 0 lt 0

###############################################################################
reset

load "myStyleEps.gnu"

print " align sum4p"

#
set output "../eps/rp/align_sum4p.eps"

set xrange [-0.50:0.50]
set yrange [-0.15:0.15]

set xlabel "{/Symbol S}_4 p_x [GeV]"
set ylabel "{/Symbol S}_4 p_y [GeV]"

#
set term post eps enh color dashed dl 2 "Helvetica" 30 size 2*5in,3*4.5in
set multiplot layout 2,3 margins 0.05,0.95,0.05,0.95

f(x,y) = a*exp(- (x-x0)**2/2/sx**2 - (y-y0)**2/2/sy**2)

a = 1e+2
x0 = 0.01; sx = 0.05
y0 = 0.01; sy = 0.05

x0_err = 0
y0_err = 0

set print "../out/rp/sum4p.par"

n = 2
do for [topo in "TB BT TT BB"] {
fit f(x,y) sprintf("../out/rp/%dpart/sum4p_%s.his",n,topo) \
    u 1:2:3:4 zerror via a, x0,sx,y0,sy

 print topo," ",n," ",a, x0,x0_err, y0,y0_err, sx,sy

 plot sprintf("../out/rp/%dpart/sum4p_%s.his",n,topo) w image
}

n = 0
do for [topo in "TB BT"] {
fit f(x,y) sprintf("../out/rp/%dpart/sum4p_%s.his",n,topo) \
    u 1:2:3:4 zerror via a, x0,sx,y0,sy

 print topo," ",n," ",a, x0,x0_err, y0,y0_err, sx,sy

 plot sprintf("../out/rp/%dpart/sum4p_%s.his",n,topo) w image
}

unset multiplot

set print
set auto y

#
set xtics ("TB" 0, "BT" 1, "TT" 2, "BB" 3)

set xlabel "Roman pot trigger topology"

c = 1e+3

set output "../eps/rp/sum4p_mean.eps"

set term post eps enh color dashed dl 2 "Helvetica" 30 size 2*5in,4.5in
set multiplot layout 1,2 margins 0.08,0.98,0.15,0.98

set size square

set key bottom right

loc(topo) = (topo eq "TB" ? 0 : \
            (topo eq "BT" ? 1 : \
            (topo eq "TT" ? 2 : 3)))

set yrange [-25:25]
set ylabel "Mean {/Symbol S}p_x [MeV]"

plot [-1:4] \
 "<awk '{if($2==0) print}' ../out/rp/sum4p.par" \
   u (loc(strcol(1))):(c*$4):(c*$5) t "0 part" w e ls 1, \
 "<awk '{if($2==2) print}' ../out/rp/sum4p.par" \
   u (loc(strcol(1))):(c*$4):(c*$5) t "2 part" w e ls 3, \
 0 lt 0

set ylabel "Mean {/Symbol S}p_y [MeV]"
plot [-1:4] \
 "<awk '{if($2==0) print}' ../out/rp/sum4p.par" \
   u (loc(strcol(1))):(c*$6):(c*$7) t "0 part" w e ls 1, \
 "<awk '{if($2==2) print}' ../out/rp/sum4p.par" \
   u (loc(strcol(1))):(c*$6):(c*$7) t "2 part" w e ls 3, \
 0 lt 0

unset multiplot

