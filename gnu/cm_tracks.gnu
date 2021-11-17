load "myStyleEps.gnu"

#
unset colorbox
set xlabel "{/Symbol h}"
set ylabel "p_T [GeV]"

set output "../eps/track/track_etaPtPos.eps"
plot [][0:1.2]   "../out/track/etaPtPos.his" u 1:2:3 w image

set output "../eps/track/track_etaPtNeg.eps"
plot [][0:1.2]   "../out/track/etaPtNeg.his" u 1:2:3 w image

#
set xlabel "{/Symbol h}"
set ylabel "{/Symbol f}"

set output "../eps/track/track_etaPhiPos.eps"
plot [][-pi:pi]  "../out/track/etaPhiPos.his" u 1:2:3 w image

set output "../eps/track/track_etaPhiNeg.eps"
plot [][-pi:pi]  "../out/track/etaPhiNeg.his" u 1:2:3 w image

set ylabel
#

set offsets graph 0, graph 0, graph 0.05, graph 0

set style data his
set ylabel "Density" offset 1.5,0

set yrange [0:]

f(x) = a*exp(-((x-m)/s)**2/2)
a = 1e-1; s = 5; m = 1e-3

#
set output "../eps/track/track_z.eps"
set xlabel "z [cm]"
fit f(x) "../out/track/trk_z.his" u 1:2 via a,m,s
set arrow 1 from -15,graph 0.2 to -15,graph 0.1 lt 7
set arrow 2 from  15,graph 0.2 to  15,graph 0.1 lt 7

set label 1 sprintf("m_z = %.2f cm",m)           at graph 0.95,0.8 right
set label 2 sprintf("{/Symbol s}_z = %.2f cm",s) at graph 0.95,0.7 right

plot [-20:20] f(x) ls 33, "../out/track/trk_z.his" u 1:2 ls 1

unset label 1
unset label 2

# Cauchy
f(x) = a/(1 + (x/(G/2))**2)
a = 10; G = 0.05

#
set output "../eps/track/track_dz.eps"

fit [-0.5:0.5] f(x) "../out/track/trk_dz.his" u 1:2 via a,G
G = abs(G)

set xlabel "z_4 - z_3 [cm]"
set arrow 1 from -10,graph 0.2 to -10,graph 0.1 lt 7
set arrow 2 from  10,graph 0.2 to  10,graph 0.1 lt 7
set label 1 sprintf("{/Symbol G} = %.3f cm",G) at graph 0.9,0.7 right
plot [-2:2] f(x) ls 33, "../out/track/trk_dz.his" u 1:2 ls 1
unset label 1

#
set output "../eps/track/track_y.eps"
set xlabel "y"
set ytics 0,0.1
plot [-3:3] \
  "../out/track/trk_y_pi.his" u 1:2 t "{/Symbol p}" ls 1, \
  "../out/track/trk_y_ka.his" u 1:2 t "K" ls 2, \
  "../out/track/trk_y_pr.his" u 1:2 t "p" ls 3
set ytics auto

set output "../eps/track/reso_y.eps"
set xlabel "y"
set ytics 0,0.1
plot [-3:3] \
  "../out/track/reso_y_pi.his" u 1:2 t "{/Symbol p}^+{/Symbol p}^{/Symbol -}" ls 1, \
  "../out/track/reso_y_ka.his" u 1:2 t "K^+K^{/Symbol -}" ls 2, \
  "../out/track/reso_y_pr.his" u 1:2 t "p@^{/Symbol -}p" ls 3
set ytics auto

#
set output "../eps/track/track_t.eps"

fit [-0.2:0.2] f(x) "../out/track/trk_t.his" u 1:2 via a,G
G = abs(G)

set xlabel "Transverse impact parameter [cm]"
set arrow 1 from -2,graph 0.2 to -2,graph 0.1 lt 7
set arrow 2 from  2,graph 0.2 to  2,graph 0.1 lt 7
set label 1 sprintf("{/Symbol G} = %.3f cm",G) at graph 0.9,0.7 right
plot [-0.5:0.5] f(x) ls 33, "../out/track/trk_t.his" u 1:2 ls 1
unset label 1

unset arrow

set style line 2 lt 2 lw 3 lc @green pt  8 dt 2
set style line 3 lt 3 lw 3 lc 3      pt  4 dt 3
set style line 4 lt 4 lw 3 lc 4      pt 12 dt 4
set style line 5 lt 5 lw 3 lc 7      pt 10 dt 5


#
#set output "../eps/track/track_sp.eps"
#set xlabel "|p@^{/ZapfDingbats\333}_3 + p@^{/ZapfDingbats\333}_4|/m"
#set arrow from 0.2,graph 0.2 to 0.2,graph 0.0
#set label 25 "cut" at 0.22,graph 0.22 left
#
#plot [0:2] \
#  "../out/track/trk_sp_pi.his" u 1:2 t "{/Symbol p}" ls 1, \
#  "../out/track/trk_sp_ka.his" u 1:2 t "K" ls 2, \
#  "../out/track/trk_sp_pr.his" u 1:2 t "p" ls 3
#unset label 25

#
f(k,x) = (x**(k/2. - 1) * exp(-x/2.)) / (2**(k/2.) * gamma(k/2.))

set xlabel "{/Symbol c}^2"

set xrange [0:35]
set yrange [0:0.2]

#
parts = "pi ka pr"
names = "'{/Symbol p}' K p"

do for [i=1:words(parts)] {
 p = word(parts,i)
 n = word(names,i)

 set label 5 n at graph 0.075,0.9 front left font "Helvetica,40"

 #
 set output sprintf("../eps/track/track_chi2a_%s.eps",p)

 plot \
 f( 2,x) t "ndf = 2"  ls 1, \
 f( 6,x) t "ndf = 6"  ls 2, \
 f(10,x) t "ndf = 10" ls 3, \
 f(14,x) t "ndf = 14" ls 4, \
 f(18,x) t "ndf = 18" ls 5, \
 sprintf("../out/track/trk_chi2_%s.his",p) \
    u ($1== 2 ? $2 : 1/0):3 ls 1 dt 1, \
 "" u ($1== 6 ? $2 : 1/0):3 ls 2 dt 1, \
 "" u ($1==10 ? $2 : 1/0):3 ls 3 dt 1, \
 "" u ($1==14 ? $2 : 1/0):3 ls 4 dt 1, \
 "" u ($1==18 ? $2 : 1/0):3 ls 5 dt 1

 #
 set output sprintf("../eps/track/track_chi2b_%s.eps",p)

 plot \
 f( 4,x) t "ndf = 4"  ls 1, \
 f( 8,x) t "ndf = 8"  ls 2, \
 f(12,x) t "ndf = 12" ls 3, \
 f(16,x) t "ndf = 16" ls 4, \
 f(20,x) t "ndf = 20" ls 5, \
 sprintf("../out/track/trk_chi2_%s.his",p) \
    u ($1== 4 ? $2 : 1/0):3 ls 1 dt 1, \
 "" u ($1== 8 ? $2 : 1/0):3 ls 2 dt 1, \
 "" u ($1==12 ? $2 : 1/0):3 ls 3 dt 1, \
 "" u ($1==16 ? $2 : 1/0):3 ls 4 dt 1, \
 "" u ($1==20 ? $2 : 1/0):3 ls 5 dt 1
}

