load "myStyleEps.gnu"

set term post eps enh color dashed dl 2 "Helvetica" 30 size 5in,4.5in

reg = "pi-low"

#
set output sprintf("../eps/res/rp_phi_c_%s.eps",reg)
set xlabel "-(t_1+t_2) [GeV^2]"
set ylabel "c [{/Symbol @{\140\140}\326}nb/GeV]"

k = 0.

# set log y ; set format y "10^{%T}"

f(x) = a*exp(-b*x)
fit f(x) sprintf("../out/res/rp_phi_res_%s.dat",reg) \
     u ($1**2+$2**2):($7/(4*$1*$2)**k/($1!=$2 ? 2 : sqrt(2))):\
                     ($8/(4*$1*$2)**k/($1!=$2 ? 2 : sqrt(2))) zerror via a,b

set label 8 \
 sprintf("d = %.2f {/Symbol\261} %.2f GeV^{-2}",\
 b,b_err) at graph 0.3,0.8

set label 9 \
 sprintf("c_0 = %.2f {/Symbol\261} %.2f {/Symbol @{\140\140}\326}nb/GeV",\
 a,a_err) at graph 0.3,0.9
 
plot sprintf("../out/res/rp_phi_res_%s.dat",reg) \
  u ($1**2+$2**2):($8<abs($7) ? \
       ($7/(4*$1*$2)**k/($1!=$2 ? 2 : sqrt(2))) : 1/0):\
       ($8/(4*$1*$2)**k/($1!=$2 ? 2 : sqrt(2))) w e ls 1 lc 0,\
     f(x) ls 1

unset label 8
unset label 9

# unset log y ; set format y

#
set output sprintf("../eps/res/rp_phi_A_%s.eps",reg)

set xlabel "-(t_1+t_2) [GeV^2]"
set log y ; set format y "10^{%T}"
set ylabel "A/(4 p_{1,T} p_{2,T}) [{/Symbol @{\140\140}\326}nb/GeV^3]" offset 0.5,0

f(x) = a*exp(-b*x)
fit [0.:] f(x) sprintf("../out/res/rp_phi_res_%s.dat",reg) \
     u ($1**2+$2**2):($3 / (4*$1*$2) / ($1!=$2 ? 2 : sqrt(2))):\
                     ($4 / (4*$1*$2) / ($1!=$2 ? 2 : sqrt(2))) zerror via a,b

set label 8 \
 sprintf("b = %.2f {/Symbol\261} %.2f GeV^{-2}",\
 b,b_err) at graph 0.1,0.1

set label 9 \
 sprintf("A_0 = %.2f {/Symbol\261} %.2f {/Symbol @{\140\140}\326}nb/GeV^3",\
 a,a_err) at graph 0.1,0.2

lc(x,y) = (x < 0.7 && y < 0.7 ? 0 : 1)
lc(x,y) = (abs(x-y) < 0.3 ? 0 : 1)

plot sprintf("../out/res/rp_phi_res_%s.dat",reg) \
     u (abs($4/$3) < 1 ? $1**2+$2**2 : 1/0):(\
         $3 / (4*$1*$2) / ($1!=$2 ? 2 : sqrt(2))):\
        ($4 / (4*$1*$2) / ($1!=$2 ? 2 : sqrt(2))):(0) w e ls 1 lc var, \
     f(x) ls 1

#
unset label 8
unset label 9
unset log y ; set format y

#
set output sprintf("../eps/res/rp_phi_R_%s.eps",reg)
set ylabel "R" offset 1.5,0

lc(x,y) = (x < 0.70 && y < 0.70 ? 0 : \
          (y < 0.64 ? 5 : \
          (y < 0.68 ? 4 : \
          (y < 0.72 ? 3 : \
          (y < 0.76 ? 2 : 1)))))

lc(x,y) = (abs(x-y) < 0.35 ? 0 : 3)

#set label "0.76 < p_{2,T} < 0.80 GeV" at graph 0.95,0.40 right tc lt 1
#set label "0.72 < p_{2,T} < 0.76 GeV" at graph 0.95,0.34 right tc lt 2
#set label "0.68 < p_{2,T} < 0.72 GeV" at graph 0.95,0.28 right tc lt 3
#set label "0.64 < p_{2,T} < 0.68 GeV" at graph 0.95,0.22 right tc lt 4
#set label "0.60 < p_{2,T} < 0.64 GeV" at graph 0.95,0.16 right tc lt 5
#set label "       p_{2,T} < 0.60 GeV" at graph 0.95,0.10 right tc lt 0
 
plot [][-2.2:1] sprintf("../out/res/rp_phi_res_%s.dat",reg) \
     u (abs($6/$5) < 0.5 ? $1**2+$2**2 : 1/0):5:6:(lc($1,$2)) w e ls 1 lc var, 0 lt 0

