all :   .rp_double .rp_occupancy .rp_patterns .rp_planeShifts \
	.rp_stripEffic .rp_stripEffic_vs_run .rp_groupEffic \
	.rp_location .rp_location_nf .rp_vtx .rp_proton .rp_weight \
	.rp_resolutions .rp_simulation .rp_veto \
	.cm_combined .cm_massResol .cm_pid .cm_relsig .cm_pid_demo \
	.cm_tracks .cm_effic .cm_ptResol \
	.lumi \
	.selection \
	.res_vir .res_mas .res_quality .res_phi_sum \
	.tunes .tunes_fac

##############################################################################
.tunes : tunes.gnu tunes.awk $(wildcard ../tune/tunes_*/results.txt) soft.dat
	gnuplot $<
	@ touch $@

.tunes_fac : tunes_fac.gnu $(wildcard ../tune/tunes_*/results.txt)
	gnuplot $<
	@ touch $@

##############################################################################
.phis : phis.gnu $(wildcard ../res/phis_*.his)
	gnuplot $<
	@ touch $@

.res_vir : res_vir.gnu $(wildcard ../out/res/hat_*.his) \
                       $(wildcard ../out/res/sha_*.his) \
         $(wildcard ../tune/tunes_*/*.yoda) \
         $(wildcard ../tune/best/*_histos.yoda)
	gnuplot $<
	@ touch $@

.res_mas : res_mas.gnu $(wildcard ../out/res/mas_*.his) \
                       $(wildcard ../out/res/mph_*.his) \
         $(wildcard ../tune/tunes_*/*.yoda) \
         $(wildcard ../tune/best/*_histos.yoda)
	gnuplot $<
	@ touch $@

.res_quality : res_quality.gnu res_quality.awk \
               $(wildcard ../out/res/rp_phi_pi-low_*.his)
	gnuplot $<
	@ touch $@

.res_phi_sum : res_phi_sum.gnu res_phi_R.gnu .res_phi
	gnuplot res_phi_sum.gnu
	gnuplot res_phi_R.gnu
	@ touch $@

.res_phi : res_phi.gnu $(wildcard ../out/res/rp_phi_*.his) \
                       $(wildcard ../tune/tunes_*/*.yoda) \
                       $(wildcard ../tune/best/*_histos.yoda)
	gnuplot $<
	@ touch $@

###############################################################################
.lumi : lumi.gnu lumi.awk ../out/lumi/lumi.out ../out/lumi/lumi.dat
	gnuplot $<
	@ touch $@

###############################################################################
.selection : selection.gnu myStyleEps.gnu \
		$(wildcard ../out/sel/?part/*.his) \
                $(wildcard ../out/rp/?part/sum?p_*.his)
	gnuplot $<
	@ touch $@

###############################################################################
.cm_combined : cm_combined.gnu cm_combined.awk \
               $(wildcard ../out/track/cmEff/*.his.gz)
	gnuplot $<
	@ touch $@

###############################################################################
.cm_pid_demo : cm_pid_demo.gnu $(wildcard ../out/track/pidDemo_*.out)
	gnuplot $<
	@ touch $@

.cm_relsig : cm_relsig.gnu $(wildcard ../out/track/relSigma_*.his)
	gnuplot $<
	@ touch $@

.cm_pid : cm_pid.gnu $(wildcard ../out/track/eloss*.his) \
	  _cm_pid.gnu cm_epsilon.gnu
	gnuplot $<
	@ touch $@

###############################################################################
.rp_resolutions : rp_resolutions.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/0part/sum2px_vs_pred_*.his) \
		$(wildcard ../out/rp/0part/sum2py_*.his) \
                $(wildcard ../out/rp/?part/sum4p_*.his)
	gnuplot $<
	@ touch $@

.rp_exclusivity : rp_exclusivity.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/sum*p_*.his) \
		$(wildcard ../out/rp/*/sump*_*.his)
	gnuplot $<
	@ touch $@

.rp_proton : rp_proton.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/prot/*.his)
	gnuplot $<
	@ touch $@

.rp_weight : rp_weight.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/weight_*.his)
	gnuplot $<
	@ touch $@

.rp_location_nf : rp_location_nf.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/hits/nf_*.his)
	gnuplot $<
	@ touch $@

.rp_location : rp_location.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/hits/loc_*.his)
	gnuplot $<
	@ touch $@

.rp_vtx : rp_vtx.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/vtx/*.his)
	gnuplot $<
	@ touch $@

.rp_simulation : rp_simulation.gnu \
         $(wildcard ../out/rp/eff_*.his) ../out/rp/coverage.his
	gnuplot $<
	@ touch $@

.rp_veto : rp_veto.gnu $(wildcard ../out/rp/2part/veto_py.his)
	gnuplot $<
	@ touch $@

.rp_double : rp_double.gnu rp_double.awk myStyleEps.gnu \
		$(wildcard ../out/rp/*/patterns.out.gz)
	gnuplot $<
	@ touch $@

.rp_occupancy : rp_occupancy.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/occupancy.out.gz)
	gnuplot $<
	@ touch $@

.rp_patterns : rp_patterns.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/patterns.fit.gz)
	gnuplot $<
	@ touch $@

.rp_planeShifts : rp_planeShifts.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/planeShifts*.dat)
	gnuplot $<
	@ touch $@

.rp_stripEffic : rp_stripEffic.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/stripEffic/*.dat.gz)
	gnuplot $<
	@ touch $@

.rp_stripEffic_vs_run : rp_stripEffic_vs_run.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/stripEffic/*.dat.gz)
	gnuplot $<
	@ touch $@

.rp_groupEffic : rp_groupEffic.gnu myStyleEps.gnu \
		$(wildcard ../out/rp/*/groupEffic/*.his)
	gnuplot $<
	@ touch $@

###############################################################################
list_part = pip pim kap kam prp prm

.cm_massResol : cm_massResol.gnu $(wildcard ../out/track/massResol/*.dat)
	gnuplot $<
	@ touch $@

.cm_ptResol : cm_ptResol.gnu .cm_fit_ptResol
	gnuplot $<
	@ touch $@

.cm_fit_ptResol   : $(foreach i,$(list_part),.cm_fit_ptResol_$(i))
	@ touch $@

# output ../out/track/ptResol_*.dat
.cm_fit_ptResol_% : cm_fitPtResol.gnu ../out/track/ptResol_%.his
	gnuplot -e "p = '$*'" $<
	@ touch $@

.cm_tracks : cm_tracks.gnu $(wildcard ../out/track/trk_*.his) \
                           $(wildcard ../out/track/reso_*.his)
	gnuplot $<
	@ touch $@

# based on simulation
.cm_effic : cm_effic.gnu $(wildcard ../out/track/sim_eff*_*.his)
	gnuplot $<
	@ touch $@
