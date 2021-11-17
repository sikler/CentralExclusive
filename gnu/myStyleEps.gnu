set key top right Left reverse samplen 3 width +0 box noauto 
set key nobox
set bar small

# set timestamp "NEW" font "Helvetica,50" textcolor lt 3
# set timestamp "%d/%m %H:%M" bottom font "Helvetica" tc rgb "dark-gray"

set fit quiet

# Margin
set lmargin at screen 0.175
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.95

# set offset 0,0, graph 0.1,0

# Terminal
set term post eps enh color dashed dl 2 "Helvetica" 30 size 5in,4.5in

set pointsize 1.5
set tics scale 2,1 front

set mxtics 5
set mytics 5

set ylabel offset 1,0

set macros
set fit logfile "/dev/null" errorvariables

set samples 1000

# Styles, colors
green = "rgb \"dark-green\""

set style line 1 lt 1 lw 3 lc 1      pt  6 dt 1  ps 1
set style line 2 lt 2 lw 3 lc @green pt  8 dt 2  ps 1
set style line 3 lt 3 lw 3 lc 3      pt  4 dt 3  ps 1
set style line 4 lt 4 lw 3 lc 4      pt 12 dt 4  ps 1
set style line 5 lt 5 lw 3 lc rgb "light-blue"      pt 10 dt 5  ps 1
set style line 6 lt 6 lw 3 lc 5      pt 14 dt 5  ps 1

#set style line 2 lt 1 lw 3 lc rgb "light-green"   pt  8 dt 2
#set style line 3 lt 1 lw 3 lc rgb "royalblue"     pt  4 dt 3
#set style line 4 lt 1 lw 3 lc rgb "light-magenta" pt 12 dt 4
#set style line 5 lt 1 lw 3 lc 7                   pt 10 dt 5

set style line 33 lt 3 lw 5 lc 3     pt  4 dt 1
set style line 55 lt 1 lw 5 lc 9     pt  4 dt 1

# rgb
gray_lambda = 0.015

color1(gray) = sqrt(gray)     + exp(-gray/gray_lambda)
color2(gray) = gray**3        + exp(-gray/gray_lambda)
color3(gray) = sin(2*pi*gray) + exp(-gray/gray_lambda)

set palette model RGB functions color1(gray), color2(gray), color3(gray)

set colorbox front
set colorbox vertical user origin 0.825,0.2 size 0.03,0.3 front

#
pip = "\"{/Symbol p}^+\""
pim = "\"{/Symbol p}^{/Symbol -}\""
kap = "\"K^+\""
kam = "\"K^{/Symbol -}\""
prp = "\"p\""
prm = "\"@^{/Symbol -}p\""

#
pipi = "\"{/Symbol p}^+{/Symbol p}^{/Symbol -}\""
kaka = "\"K^+K^{/Symbol -}\""
prpr = "\"p@^{/Symbol -}p\""

# index FIXME
dime_ = 0
supc_ = 1
gran_ = 2
mine_ = 3

# title
dime = "\"Dime\""
supc = "\"SuperChic3\""
gran = "\"Graniitti\""
mine = "\"Our param.\""
data = "\"Data\""

# ls
_dime = 2
_supc = 4
_gran = 3
_mine = 1
_data =-1

# masses
mpi = 0.139570
mka = 0.493677
mpr = 0.938272
mel = 0.510999e-3

#
lab(i) = sprintf("%d%s%s%s", (i&8 ? 2 : 1), \
                             (i&4 ? "f" : "n"), \
                             (i&2 ? "B" : "T"), \
                             (i&1 ? "v" : "u"))
