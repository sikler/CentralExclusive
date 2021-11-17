load "myStyleEps.gnu"

set xrange [0:pi]
set yrange [0:pi]
set zrange [0:]
set cbrange [0:]

set style data image

set xlabel "@^{\\176}{/Symbol f}"
set ylabel "{/Symbol f}"

set size square

unset colorbox

#
f(x) = pi * (1 - (1 - x/pi)**n)**(1/n)

#
g(x) = pi * (1 - (1 - (x/pi)**n)**(1/n))

set label 1 at graph 0.9,0.1 right front

set output "../eps/phis_pi.eps"
n = 1.11
set label 1 sprintf("n = %.2f",n)
plot  "../res/phis_pi.his" u 1:2:3, f(x) ls 2 lc 2, x ls 4

set output "../eps/phis_ka.eps"
n = 1.05
set label 1 sprintf("n = %.2f",n)
plot  "../res/phis_ka.his" u 1:2:3, f(x) ls 2 lc 2, x ls 4

set output "../eps/phis_pr.eps"
n = 1.03
set label 1 sprintf("n = %.2f",n)
plot  "../res/phis_pr.his" u 1:2:3, f(x) ls 2 lc 2, x ls 4

