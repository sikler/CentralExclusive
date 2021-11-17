load "myStyleEps.gnu"

nChi = 3.4 # FIXME
nTop = 5.2 # FIXME

c = 1e-3

#
set xlabel "x@_1* [cm]"
set ylabel "x@_2* [cm]" offset 2,0
set size square
set xrange [-0.1:0.1]
set yrange [-0.1:0.1]

set label 6 at graph 0.05,0.9 left front
set label 7 at graph 0.05,0.8 left front

set style data image
#set cblabel "counts [10^3]" rotate parallel
#set colorbox user origin 0.7,0.2
unset colorbox

do for [n in "0 2"] {

if(n eq "0") {
  set label 6 "pp"
} else {
  set label 6 "p(h^+h^-)p"
}

set label 7 "elastic"
set output sprintf("../eps/sel/rp_x12_elas_%spart.eps",n)
plot       sprintf("../out/sel/%spart/rp_x12_elas.his",n) u 1:2:(c*$3), \
 x-0.01 lt 2, x+0.01 lt 2

set label 7 "exclusive"
set output sprintf("../eps/sel/rp_x12_sign_%spart.eps",n)
plot       sprintf("../out/sel/%spart/rp_x12_sign.his",n) u 1:2:(c*$3), \
 x-0.01 lt 2, x+0.01 lt 2

set label 7 "sideband"
set output sprintf("../eps/sel/rp_x12_side_%spart.eps",n)
plot       sprintf("../out/sel/%spart/rp_x12_side.his",n) u 1:2:(c*$3), \
 x-0.01 lt 2, x+0.01 lt 2
}

unset label 6
unset label 7

set auto x
set auto y
set size nosquare

#
set output "../eps/sel/cm_z.eps"
set xlabel "z [cm]"
set ylabel "counts [10^3]" offset -0.5,0

f(x) = a*exp(-abs((x-m)/s)**2/2) + b
a = 1; m = 0.1; s = 5; b = 10
fit f(x)   "../out/sel/2part/cm_z.his" u 1:(c*$2):(c*$3) zerror via a,m,s,b
s = abs(s)

set log y; set format y "10^{%T}"
set label 6 sprintf("m_z = %.2f cm\n{/Symbol s}_z = %.2f cm",m,s) \
    at graph 0.95,0.9 right
plot [][1e-2:] "" u 1:(c*$2) w his ls 1, f(x)-b ls 3
unset log y; set format y
print sprintf(" tracks  : z_mean = %.2f cm, z_sigma = %.2f [cm]",m,s)
unset label 6

# 
set output "../eps/sel/cm_dz.eps"
set xlabel "z_1 - z_2 [cm]"
set log y; set format y "10^{%T}"
f(x) = a*exp(-x**2/2/s**2)
a = 1
s = 5*sqrt(2)
fit f(x)    "../out/sel/2part/cm_dz.his" u 1:(c*$2):(abs($1) > 5 ? c*$3 : 1e+30) \
            zerror via a,s
plot [][1e-2:] "" u 1:(c*$2) w his ls 1 #, (abs(x) > 5 ? f(x) : 1/0) ls 3
unset log y; set format y
print sprintf(" distance: dz_sigma = %.2f [cm]",s)

#
set output "../eps/sel/cm_dzrel.eps"
f(x) = a*exp(-abs(x/s)**2/2)
a = 1; s = 1

set ylabel offset 0,0

fit [-5:5] f(x) "../out/sel/2part/cm_dzrel.his" u 1:(c*$2):(c*$3) zerror via a
set xlabel "(z_1 - z_2) / ({/Symbol s}@_1^2 + {/Symbol s}@_2^2)^{1/2}"
plot [-8:8][0:]          "" u 1:(c*$2) w his ls 1, \
     f(x) ls 3
unset log y; set format y
print sprintf(" distance: dzrel_sigma = %.2f (set)",s)

#
set log cb; set format cb "10^{%T}"

set output "../eps/sel/cm_xy.eps"
set xlabel "x_{vtx,2h} [cm]"
set ylabel "y_{vtx,2h} [cm]"
set object 1 circle at 0.1,0.1 radius 1 fs empty border lt 2 front
set size square
plot       "../out/sel/2part/cm_xy.his" w image
set size nosquare
unset object 1

set output "../eps/sel/cm_zr.eps"
set xlabel "z_{vtx,2h} [cm]"
set ylabel "r_{vtx,2h} [cm]"
plot       "../out/sel/2part/cm_zr.his" w image, 1 lt 2

unset log cb

#
unset colorbox
set xtics out
set ytics out

set output "../eps/sel/chi_24.eps"
set xlabel "{/Symbol c} [pp]" offset 0,0.5
set ylabel "{/Symbol c} [p(h^+h^-)p]" offset 0,0

set size square

set arrow from 0,0 to nChi,nChi                lt 2 nohead front
set arrow from nChi,nChi to graph 1,first nChi lt 2 nohead front
set arrow from nChi,nChi to nChi,graph 1       lt 2 nohead front

set arrow from nTop,nTop to graph 1,first nTop lt 2 dt 2 nohead front

plot       "../out/sel/2part/chi_24.his" w image #, (x < 3 ? x : 3) lt 2

set xtics in
set ytics in
unset arrow

# for out to
m = 12

pdf(x) = x*exp(-x**2/2)
bac(x) = abs(b) * x * exp(-k*abs(x)**d)

f(x) = (x < m ? abs(a) * (pdf(x) + bac(x)) : 1/0)
g(x) = (x < m ? abs(a) * (         bac(x)) : 1/0)

a = 1; b = 0.05; k = 0.1; d = 1

set arrow from nChi,graph 0 rto 0,graph 1 lt 2 nohead
set arrow from nTop,graph 0 rto 0,graph 1 lt 2 dt 2 nohead

set size nosquare
set xlabel offset 0,0
set ylabel "counts [10^3]" offset 1,0

set label 6 at graph 0.95,0.2 right

# chi_2p
set output "../eps/sel/chi_2.eps"
set xlabel "{/Symbol c} [pp]"

fit [0:m] f(x) \
           "../out/sel/2part/chi_2.his" u 1:(c*$2):(c*$3) zerror via a
fit [0:m] f(x) \
           "../out/sel/2part/chi_2.his" u 1:(c*$2):(c*$3) zerror via a,b,k

set label 6 sprintf("k = %.3f",k)
plot [0:m] "../out/sel/2part/chi_2.his" u 1:(c*$2) w his ls 1, \
           f(x) ls 3, g(x) ls 2
#          "../out/sel/2part/oth_4.his" u 1:(c*$2/10 * 1.3) w his lt 2
print sprintf(" chi_2p : k = %.3f",k)

# chi_4p
set output "../eps/sel/chi_4.eps"
set xlabel "{/Symbol c} [p(h^+h^-)p]"

fit [0:m] f(x) \
           "../out/sel/2part/chi_4.his" u 1:(c*$2):(c*$3) zerror via a
fit [0:m] f(x) \
           "../out/sel/2part/chi_4.his" u 1:(c*$2):(c*$3) zerror via a,b,k

set label 6 sprintf("k = %.3f",k)

w(x) = (x < nChi ? x/nChi : 1)

plot [0:m]  \
 "../out/sel/2part/oth_4.his" u 1:($1 < 2.5 || $1 > 4.5 ? c*$2*0.20/w($1) : 1/0) w l lt 4 lw 5, \
 "../out/sel/2part/chi_4.his" u 1:(c*$2) w his ls 1, \
 f(x) ls 3, g(x) ls 2

print sprintf(" chi_4p : k = %.3f",k)

#
set xlabel
set ylabel
set term post eps enh color dashed dl 2 "Helvetica" 25 size 2*5in,2.5*4.5in

set output "../eps/sel/chi_4_dphi.eps"

set label 7 at graph 0.95,0.375 right
set label 6 at graph 0.95,0.250 right

set multiplot layout 5,4 margins 0.10,1.00,0.05,0.95 spacing 0,0.02
set size square

unset xtics
unset ytics

unset xlabel
unset ylabel

set print "../out/sel/k_vs_dphi.par"

do for [i=0:17] {
 if(i >= 14)    { set xtics ; set xlabel "{/Symbol c} [p(h^+h^-)p]" }
 if(i % 4 == 0) { set ytics ; set ylabel "counts [10^3]" offset 1,0}

 set label 7 sprintf("%d < {/Symbolf f} < %d{/Symbol\260}",i*10,(i+1)*10)
 fit [0:m] f(x) sprintf("../out/sel/2part/chi_4_%d.his",i) \
           u 1:(c*$2):(c*$3) zerror via a
 fit [0:m] f(x) "" \
           u 1:(c*$2):(c*$3) zerror via a,b,k

 set label 6 sprintf("k = %.3f",k)

 plot [0:m] "" u 1:(c*$2) w his ls 1, f(x) ls 3, g(x) ls 2
#      sprintf("../out/sel/2part/oth_4_%d.his",i) \
#        u 1:(c*$2/5./($1 < 3.4 ? ($1/3.4) : 1)) w his ls 4

 int_b(x) = x*exp(-k*x)
 int_sign = (1 - exp(-k*nChi) * (k*nChi + 1))/k**2
 int_side = 0
 x = nChi; dx = 1e-3

 while (1) {
  x = x + dx
  int_side = int_side + int_b(x-dx/2)*dx
  if(int_side > int_sign) { break }
 }

 print i*10,(i+1)*10, k,k_err, x

 unset xtics
 unset ytics

 unset xlabel
 unset ylabel
}
unset multiplot

set size noratio

###############################################################################o
reset
load "myStyleEps.gnu"

set output "../eps/sel/k_vs_dphi.eps"

#set ytics 0,0.02

set xlabel "{/Symbol f} [{/Symbol\260}]" offset 0,0.5
set ylabel "k" offset 1.5,0

plot [][0:] "../out/sel/k_vs_dphi.par" u (($1+$2)/2):3:4 w e ls 1

###############################################################################o
reset
load "myStyleEps.gnu"

set term post eps enh color dashed dl 2 "Helvetica" 25 size 2*5in,2.5*4.5in

set key noauto
unset colorbox

set output "../eps/sel/chi_24_dphi.eps"

set label 7 at graph 0.5,0.925 center front

set multiplot layout 5,4 margins 0.10,1.00,0.05,0.95 spacing 0,0.02
set size square

unset xtics
unset ytics

unset xlabel
unset ylabel

set arrow from 0,0 to nChi,nChi                lt 2 nohead front
set arrow from nChi,nChi to graph 1,first nChi lt 2 nohead front
set arrow from nChi,nChi to nChi,graph 1       lt 2 nohead front

set arrow from nTop,nTop to graph 1,first nTop lt 2 dt 2 nohead front

do for [i=0:17] {
 if(i >= 14)    { set xtics ; set xlabel "{/Symbol c} [p(h^+h^-)p]" }
 if(i % 4 == 0) { set ytics ; set ylabel "{/Symbol c} [pp]" offset 1,0}

 set label 7 sprintf("%d < {/Symbolf f} < %d{/Symbol\260}",i*10,(i+1)*10)

 plot sprintf("../out/sel/2part/chi_24_%d.his",i) w image

 unset xtics
 unset ytics

 unset xlabel
 unset ylabel
}

unset multiplot

###############################################################################o
reset
load "myStyleEps.gnu"

unset colorbox
set style data image
set size square

set term post eps enh color dashed dl 2 "Helvetica" 40 size 2*5in,2*4.5in

###############################################################################
set label 5 at graph 0.95,0.9 right front
 
set xrange [-0.50:0.50]
set yrange [-0.15:0.15]

set ytics -0.5,0.1

set object circle at 0,0 size 0.02 fc lt 2 fs empty front

#
ncm = 2
set output sprintf("../eps/sel/sum%dp.eps",ncm)

set multiplot layout 2,2 margins 0.08,0.98,0.08,0.98 spacing 0.1

set object 5 rectangle from screen 0.48,screen 0 rto 0.5,1 \
    fillc rgb "yellow" fs solid 0.1 noborder behind

set xlabel sprintf("{/Symbol S}_%dp_x [GeV]",ncm) offset 0,0.5
set ylabel sprintf("{/Symbol S}_%dp_y [GeV]",ncm) offset 3,0

do for [topo in "TB BT"] {
do for [n in "2 0"] {

 set label 5 sprintf("%s, %sh",topo,n)

 plot sprintf("../out/rp/%spart/sum%dp_%s.his",n,ncm,topo)

 unset object 5
}}

unset multiplot

#
ncm = 4
set output sprintf("../eps/sel/sum%dp.eps",ncm)

set multiplot layout 2,2 margins 0.08,0.98,0.08,0.98 spacing 0.1 columnsfirst

set xlabel sprintf("{/Symbol S}_%dp_x [GeV]",ncm) offset 0,0.5
set ylabel sprintf("{/Symbol S}_%dp_y [GeV]",ncm) offset 3,0

do for [topo in "TB BT TT BB"] {
n = 2

 set label 5 sprintf("%s, %dh",topo,n)

 plot sprintf("../out/rp/%dpart/sum%dp_%s.his",n,ncm,topo)

}

unset multiplot

###############################################################################
set xrange [-1:1]
set yrange [-1:1]

set xtics -1,0.5
set ytics -1,0.5

set label 5 at graph 0.95,0.9 right front

#
set term post eps enh color dashed dl 2 "Helvetica" 20 size 5in,2*4.5in

set output "../eps/sump.eps"

set multiplot layout 4,2 margins 0.08,0.98,0.08,0.98 rowsfirst spacing 0.05

do for [topo in "TB BT TT BB"] {
do for [dir in "x y"] {

 set xlabel sprintf("{/Symbol S}_2p_%s [GeV]",dir)
 set ylabel sprintf("{/Symbol S}_4p_%s [GeV]",dir)

 set label 5 topo

 plot sprintf("../out/rp/2part/sump%s_%s.his",dir,topo) #, 0 lt 2 dt 2

}}

unset multiplot
unset label 5

#
reset
load "myStyleEps.gnu"

set key bottom right

set output "../eps/sel/cm_looper.eps"
set xlabel "|p@^{/ZapfDingbats\333}_3 + p@^{/ZapfDingbats\333}_4|/m"
set arrow from 0.2,graph 0.2 to 0.2,graph 0.0
set label 6 "cut" at 0.22,graph 0.22 left

plot [0:2] \
  "../out/sel/2part/cm_looper_pi.his" u 1:2 t "{/Symbol p}" w his ls 1, \
  "../out/sel/2part/cm_looper_ka.his" u 1:2 t "K" w his ls 2, \
  "../out/sel/2part/cm_looper_pr.his" u 1:2 t "p" w his ls 3
unset label 6

