load "myStyleEps.gnu"
load "cm_epsilon.gnu"

set colorbox front

set key height +0.3 samplen 3

set style data image

set xlabel "p [GeV]"
set ylabel "ln({/Symbol e}/[MeV/cm])" offset 1,0

set xrange [log(0.1):log(10)]

set xtics \
 ("0.1" log(0.1), "0.2" log(0.2), "" log(0.3), "" log(0.4), \
  "0.5" log(0.5), "" log(0.6), "" log(0.7), "" log(0.8), "" log(0.9), \
  "1" log(1), "2" log(2), "" log(3), "" log(4), "5" log(5), "10" log(10))

set yrange [0.0:3.5]

ele_ = 17; set style line ele_ lw 2 dt 4 lc @green
pio_ = 11; set style line pio_ lw 2      lc 1 
kao_ = 13; set style line kao_ lw 2 dt 3 lc 4
pro_ = 15; set style line pro_ lw 2 dt 2 lc 3

elec = "\"e\""
pion = "\"{/Symbol p}\""
kaon = "\"K\""
prot = "\"p\""

fpion(p) = epsilon(p/mpi)
fkaon(p) = epsilon(p/mka)
fprot(p) = epsilon(p/mpr)
fsigm(p) = epsilon(p/msi)

fpika(p) = (fpion(p)+fkaon(p))/2
fkapr(p) = (fkaon(p)+fprot(p))/2

set colorbox vertical user origin graph 0.55,graph 0.55 front

set cbtics offset -1,0

set label 5 at screen 0.2, screen 0.25 left front

set label 5 "signal"
call "_cm_pid.gnu" "elossAl" 1

set label 5 "unidentified"
call "_cm_pid.gnu" "elossNo" 0

set label 5 "{/Symbol p}^+{/Symbol p}^{/Symbol -}"
call "_cm_pid.gnu" "elossPi" 0

set label 5 "K^+K^{/Symbol -}"

call "_cm_pid.gnu" "elossKa" 0

set label 5 "p@^{/Symbol -}p"
call "_cm_pid.gnu" "elossPr" 0

set label 5 "sideband"
call "_cm_pid.gnu" "elossSb" 1

