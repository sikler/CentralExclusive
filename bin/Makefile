ONLXPLUS := $(shell uname -n | grep lx)

ifdef ONLXPLUS
nevery = 1
else
nevery = 1
endif

# 
nproc = $(shell nproc)


list_parts = pi ka pr
list_parts_charged = pip pim kap kam prp prm

###############################################################################
# trk -> rp -> anal
all :

###############################################################################
# Analysis (anal) :
# lumi_weights -> collect -> cm_eff / pid_demo / pt_resol
#              -> phys -> lumi_postproc

anal :
	make lumi_weights
	make collect
	make -j3 cm_eff
	make -j3 pid_demo
	make pt_resol
	make phys
	make lumi_postproc

################################################
# Luminosity - post-process
# input : ../out/lumi/takenEvents.dat (through .collect)
# output: ../out/lumi/lumi.out, pred vs meas events for ../gnu/lumi.gnu
lumi_postproc :     ../lumi/runs.dat \
	 $(wildcard ../lumi/*ns_*.txt) \
	 $(wildcard ../lumi/injectionSchemes.txt) \
	            ../lumi/lumi_exceptions.par \
	            ../pars/mc_efficiency.par
	$(eval flags = -nCmTra 2 -nCpus $(nproc))
	./analysis -lumi -postProcess $(flags)

################################################
# Do Physics
# ../out/track/ptResol_{part}.dat is needed for mass response
phys :
	$(eval flags = -nCmTra 2 -nCpus $(nproc) -nEvery $(nevery))
	./analysis -data -doPhys   $(flags)

###############################################################################
# Pt resolution
# input  : ../out/track/ptResol_{part}.his
# output : ../out/track/ptResol_{part}.dat

pt_resol :
	( cd ../gnu ; make .cm_fit_ptResol )

###############################################################################
# Pid : pid_demo
# input : ../out/track/relSigma_{part}.his
# output: ../out/track/pidDemo_{part}.his

pid_demo : $(foreach i,$(list_parts),.pid_demo_$(i))

.pid_demo_% :
	$(eval flags = -nCmTra 2)
	./analysis -pidDemo $* $(flags)
	@ touch $@

################################################
# Combined HLT/tracking/PID efficiency
cm_eff : $(foreach i,$(list_parts),.cm_eff_$(i))
	echo "compressing ../out/track/cmEff/*.his..."
	gzip -f ../out/track/cmEff/*.his

# input  : ../out/track/simTable_{part}.out.gz
# input  : ../out/track/sim_effOne_{part}{p,m}.his
# input  : ../out/track/relSigma_{part}.his
# output : ../out/track/cmEff/{part}_{p1T}_{p2T}_{topo}.his.gz
.cm_eff_% :
	$(eval flags = -nCmTra 2 -nCpus $(nproc))
	./analysis -calc -cmEff $* $(flags)
	@ touch $@

################################################
# Collect from data
collect :
	$(eval flags = -nCmTra 2 -nCpus $(nproc) -nEvery $(nevery))
	./analysis -data -rpReco   $(flags)
	$(eval flags = -nCmTra 0 -nCpus $(nproc) -nEvery $(nevery))
	./analysis -data -rpReco   $(flags)

################################################
# Luminosity - calculate weights
# output: ../out/lumi/lumi.dat     : run ls Leff time
# output: ../out/lumi/int_lumi.dat : (Lrec) Leff
# output: ../out/lumi/run_lumi.dat : (Lrec)
lumi_weights :      ../lumi/runs.dat \
	 $(wildcard ../lumi/*ns_*.txt) \
	 $(wildcard ../lumi/injectionSchemes.txt) \
	            ../lumi/lumi_exceptions.par \
	            ../lumi/lumi_DCSONLY.csv.gz
	$(eval flags = -nCmTra 2 -nCpus $(nproc))
	./analysis -lumi -calcWeights $(flags)

###############################################################################
# Roman pots (rp) :
# rp_local -> rp_veto -> rp_angEff

rp :
	make rp_local
	make rp_veto
	make rp_angEff

# calculate roman pot angular effic vs (phi,p1T,p2T), based on accep and !veto
# input  : ../out/rp/2part/veto_py.his
# input  : |py| acceptance (0.190,0.655) FIXME
# output : ../out/rp/eff_{TB,BT,TT,BB}.his
# output : ../out/rp/coverage.his
rp_angEff  :
	$(eval flags = -nCmTra 2 -nCpus $(nproc) -nEvery $(nevery))
	./analysis -calc -rpAngEff $(flags)

# extract roman pots !veto efficiency (relevant or TB,BT) using TT,BB tracks
# output : ../out/rp/2part/veto_py.his
rp_veto :
	$(eval flags = -nCmTra 2 -nCpus $(nproc) -nEvery $(nevery))
	./analysis -data -rpVeto   $(flags)

# roman pots local calibration
rp_local :
	make rp_2part
	make rp_0part
	make rp_generate

rp_2part :
	$(eval flags = -nCmTra 2 -nCpus $(nproc))
	./analysis -data -rpPats   $(flags) # collecting roman pot strip patt's
	./analysis -calc -rpShifts $(flags) # optimizing shifts
	./analysis -calc -rpPats   $(flags) # fitting roman pot hit patterns
	./analysis -data -rpEffic  $(flags) # collecting roman pot strip eff's
	./analysis -calc -rpGroup  $(flags) # calculate roman pot tracklet eff's

rp_0part :
	$(eval flags = -nCmTra 0 -nCpus $(nproc))
	./analysis -data -rpPats   $(flags) # collecting roman pot strip patt's
	./analysis -calc -rpShifts $(flags) # optimizing shifts
	./analysis -calc -rpPats   $(flags) # fitting roman pot hit patterns
	./analysis -data -rpEffic  $(flags) # collecting roman pot strip eff's
	./analysis -calc -rpGroup  $(flags) # calculate roman pot tracklet eff's

rp_generate :
	$(eval flags = -nCmTra 9 -nCpus $(nproc))
	./analysis -calc -rpGene   $(flags) # generate roman pot data to file
	#
	./analysis -data -rpPats   $(flags) # collecting roman pot strip patt's
	./analysis -calc -rpShifts $(flags) # optimizing shifts
	./analysis -calc -rpPats   $(flags) # fitting roman pot hit patterns

# Print roman pot events
# input : ../data/2print/TOTEM20/diag.dat.gz
# input : ../data/2print/TOTEM40/para.dat.gz
rp_2print :
	$(eval flags = -nCmTra 2)
	./analysis -data -rpPrint $(flags)

###############################################################################
# Tracker (trk) :
# trk_eff

trk :
	make -j3 trk_eff

# Tracking efficiency (using simulation)
# input : ../simu/simulation_*.dat.gz
# output: ../out/track/sim_*.his
# output: ../out/track/simTable_*.out.gz
# output: ../out/track/simTable_*.zeros
# output: ../out/track/ptResol_*.his

trk_eff : $(foreach i,$(list_parts_charged),.trk_eff_$(i))
	echo done

.trk_eff_% : $(wildcard ../simu/simulation_*.dat.gz)
	$(eval flags = -nCmTra 2)
	./analysis -trkEff $* $(flags)
	@ touch $@

###############################################################################
prof ::
	gprof -b -p analysis gmon.out | head -20

