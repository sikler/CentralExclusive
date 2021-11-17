load "myStyleEps.gnu" # 100 lines

set datafile missing "0"
set samples 100

#
nSect = 3
array centers[nSect]
centers[1] = "0.225 0.275 0.325 0.375"; np = words(centers[1])
centers[2] = "0.425 0.475 0.525 0.575"
centers[3] = "0.625 0.675 0.725 0.775"

allcenters = sprintf("%s  %s  %s",centers[1],centers[2],centers[3])

dpt = 0.05/2 # (0.8 - 0.2)/12/2

# relative syst uncertainties
runc = 5.4e-2

print sprintf(" systematic uncertainty = %.1f%", runc*1e+2)

#
data_lin(z) = \
  sprintf("<awk 'BEGIN {p1t=%s; p2t=%s} \
           {if($1>0.30 && $1<2.50 && $2==p1t && $3==p2t) print $1,$4,$5}' \
               ../out/res/mas_%s.his", p1t,p2t,p)

data_log(z) = \
  sprintf("<(awk 'BEGIN {p1t=%s; p2t=%s} \
      {if($1>0.30 && $1<1.60 && $2==p1t && $3==p2t && $4>0) print $1,$4,$5}' \
               ../out/res/mas_%s.his ; \
              awk 'BEGIN {p1t=%s; p2t=%s} \
      {if($1>1.60 && $1<2.50 && $2==p1t && $3==p2t && $4>0) print $1,$4,$5}' \
               ../out/res/mhi_%s.his )", p1t,p2t,p, p1t,p2t,p)

data(z) = (z == 1 ? data_lin(z) : data_log(z))

#
dime(fac) = \
 sprintf("<(yoda2flat %s ../tune/best/%s_%s_histos.yoda - | \
            awk '{if(NF>1 && ($1+$2)/2<1.8 && $3!=0 && $4!=0) print}' ; \
            yoda2flat %s ../tune/best/%s_%s_histos.yoda - | \
            awk '{if(NF>1 && ($1+$2)/2>1.8 && $3!=0 && $4!=0) print}')",\
            m,re,fac, m_,re,fac)

#
setXtics = \
 "set xtics \
  ('0' 0,'0.5' 0.5,'1' 1,'1.5' 1.5,'2' 2,'2.5' 2.5,'3' 3,'3.5' 3.5,'4' 4); \
  set xtics add ('' 0.1 2, '' 0.2 2, '' 0.3 2, '' 0.4 2); \
  set xtics add ('' 0.6 2, '' 0.7 2, '' 0.8 2, '' 0.9 2); \
  set xtics add ('' 1.1 2, '' 1.2 2, '' 1.3 2, '' 1.4 2); \
  set xtics add ('' 1.6 2, '' 1.7 2, '' 1.8 2, '' 1.9 2); \
  set xtics add ('' 2.1 2, '' 2.2 2, '' 2.3 2, '' 2.4 2); \
  set xtics add ('' 2.6 2, '' 2.7 2, '' 2.8 2, '' 2.9 2); \
  set xtics add ('' 3.1 2, '' 3.2 2, '' 3.3 2, '' 3.4 2); \
  set xtics add ('' 3.6 2, '' 3.7 2, '' 3.8 2, '' 3.9 2)"

unsetXtics = \
 "set xtics ('' 0, '' 0.5, '' 1, '' 1.5, '' 2, '' 2.5, '' 3); \
  set xtics add ('' 0.1 2, '' 0.2 2, '' 0.3 2, '' 0.4 2); \
  set xtics add ('' 0.6 2, '' 0.7 2, '' 0.8 2, '' 0.9 2); \
  set xtics add ('' 1.1 2, '' 1.2 2, '' 1.3 2, '' 1.4 2); \
  set xtics add ('' 1.6 2, '' 1.7 2, '' 1.8 2, '' 1.9 2); \
  set xtics add ('' 2.1 2, '' 2.2 2, '' 2.3 2, '' 2.4 2); \
  set xtics add ('' 2.6 2, '' 2.7 2, '' 2.8 2, '' 2.9 2); \
  set xtics add ('' 3.1 2, '' 3.2 2, '' 3.3 2, '' 3.4 2); \
  set xtics add ('' 3.6 2, '' 3.7 2, '' 3.8 2, '' 3.9 2)"

#
sysunc = \
 'data(z)  u 1:2:(0.02/2):($2*runc) w boxxye lc 0 fs solid 0.2 noborder'

#
result = \
 'data(z) u 1:2:3 lw 3 lt 5 lc 0 pt 4 ps 2, \
       "" u 1:2 w l lt 1 lc 0'

#
# w = 1e+4
lw = 10

dexp(lam,soft) = sprintf("{/Symbol L}=%.1f, s%d",lam,soft)

#
set style data e

###############################################################################
set term post eps enh color dashed dl 0.5 "Helvetica" 50 size np*5in,np*4.5in

aver(p1t,p2t) = exp(-4 * (p1t**2 + p2t**2))

do for [z=1:2] {
# do for [p in "pi ka"] { FIXME VERY
do for [p in "pi"] {
 print " res_mas ",p,(z == 1 ? "" : " log")

 if(z == 1) {
  if(p eq "pi") { set xrange [0.0:1.6] ; cy = 5.0 } 
  if(p eq "ka") { set xrange [0.5:2.1] ; cy = 0.5 }
  if(p eq "pr") { set xrange [1.5:3.5] ; cy = 0.5 }
 } else {
  if(p eq "pi") { set xrange [0.0:2.5] ; cy = 5.0 } 
  if(p eq "ka") { set xrange [0.5:2.5] ; cy = 0.5 }
  if(p eq "pr") { set xrange [1.5:3.5] ; cy = 0.5 }
 }

 if(p eq "pi") { re = "pipm" }
 if(p eq "ka") { re = "kpkm" }

 do for [i=1:nSect] {
 do for [j=i:nSect] {
 if(z == 1) {
   set output sprintf("../eps/res/mas_%s_%d_%d.eps",    p,i,j)
  } else {
   set output sprintf("../eps/res/mas_%s_log_%d_%d.eps",p,i,j)
   set log y #; set format y "10^{%T}"
  }

 if(i==3 && j==3) {
   unset border; unset tics
   plot [][1:2] 0
   set border; set tics
 } else {
 set multiplot layout np,np upwards \
      margins 0.06,1.00,0.06,0.99 spacing 0.02,0.01
 set size square

 x = 0
 do for [p1t in centers[i]] {
  y = 0
  do for [p2t in centers[j]] {
   if(p1t <= p2t) {

    if(z == 1) {
     set yrange [0:]
    } else {
     set yrange [cy*5e-1*aver(p1t,p2t) : \
                 cy*1e+2*aver(p1t,p2t)]
    }

    @unsetXtics; set xlabel; set ylabel
    if(x == 0) {
      @setXtics
      set xlabel sprintf("%.2f < p_{1,T} < %.2f GeV",p2t-dpt,p2t+dpt) \
          offset 0,0.3}
    if((i==j && x==y) || (i!=j && y==0)) {
      set ylabel sprintf("%.2f < p_{2,T} < %.2f GeV",p1t-dpt,p1t+dpt) \
          offset 1.8,0 }

    if(2*x == np-1 && 2*y == np-1) {
      set label 50 "m [GeV]" at graph 0.5,0.925 center front
      set label 51 "d^3{/Symbol s}/dmdp_{1,T}dp_{2,T} [nb/GeV^5]" \
                             at graph 0.075,0.5 rotate by 90 center front
    }

    if(z == 1) { unset key }
    if(z == 2) { unset key }

    if(x == 3 && y == 3) {
      set label 12 "" at graph 0.9,0.9 font "Helvetica,60" right
      if(p eq "pi") { set label 12 @pipi }
      if(p eq "ka") { set label 12 @kaka }
      if(p eq "pr") { set label 12 @prpr }
    }

     m  = sprintf("-m dndmas-%d-%d-his",4*(i-1)+x,4*(j-1)+y)
     m_ = sprintf("-m dndmhi-%d-%d-his",4*(i-1)+x,4*(j-1)+y)

     if((4*(i-1)+x <= 1 && 4*(j-1)+y <= 1)) {
       plot @result
     } else {

       if(p eq "pi") {
         c = 20e-3 # FIXME
       } else {
         c = 20e-3 * 2./3 # reduce sigo for kaons by 1./3
       }

       w(x) = 1e+6/x**2

       plot \
         dime("orear") u (($1+$2)/2):(c*($3-$4)):(c*($3+$4)) \
           w filledc fs solid 0.75 noborder t "Dime Orear" lt 2, \
         dime("exp")   u (($1+$2)/2):(c*($3-$4)):(c*($3+$4)) \
           w filledc fs solid 0.75 noborder t "Dime exp"   lt 3, \
         dime("power") u (($1+$2)/2):(c*($3-$4)):(c*($3+$4)) \
           w filledc fs solid 0.75 noborder t "Dime power" lt 1, \
         @result
     }
   
     set auto y 
     unset label 50
     unset label 51

     unset label 12
   } else { set multiplot next }
   y = y + 1
  }
  x = x + 1
 }
 unset multiplot
 }}}
}
}

