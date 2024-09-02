load glm
load dic
model in "model.txt"
data in "data.txt"
compile, nchains(1)
parameters in "inits2.txt", chain(1)
initialize
adapt 500
update 5000
monitor deviance, thin(5)
monitor grwth_Transect_randomeffect, thin(5)
monitor surv_Transect_randomeffect, thin(5)
monitor reproyesno_Transect_randomeffect, thin(5)
monitor repro_Transect_randomeffect, thin(5)
monitor grwth_intercept, thin(5)
monitor grwth_RosCoef, thin(5)
monitor grwth_TempFallCoef, thin(5)
monitor grwth_TempSummerCoef, thin(5)
monitor grwth_TempWinterCoef, thin(5)
monitor grwth_PptFallCoef, thin(5)
monitor grwth_PptSummerCoef, thin(5)
monitor grwth_PptWinterCoef, thin(5)
monitor grwthvar_intercept, thin(5)
monitor grwthvar_RosCoef, thin(5)
monitor surv_intercept, thin(5)
monitor surv_RosCoef, thin(5)
monitor surv_PptWinterCoef, thin(5)
monitor surv_TempFallCoef, thin(5)
monitor surv_TempSummerCoef, thin(5)
monitor surv_TempWinterCoef, thin(5)
monitor reproyesno_intercept, thin(5)
monitor reproyesno_RosCoef, thin(5)
monitor reproyesno_PptFallCoef, thin(5)
monitor reproyesno_PptSummerCoef, thin(5)
monitor reproyesno_TempFallCoef, thin(5)
monitor reproyesno_TempSummerCoef, thin(5)
monitor reproyesno_TempWinterCoef, thin(5)
monitor repro_intercept, thin(5)
monitor repro_RosCoef, thin(5)
monitor repro_PptFallCoef, thin(5)
monitor repro_PptSummerCoef, thin(5)
monitor repro_TempWinterCoef, thin(5)
monitor repro_TempFallCoef, thin(5)
monitor repro_TempSummerCoef, thin(5)
monitor newplt_intercept, thin(5)
monitor r.infls, thin(5)
monitor r.newplts, thin(5)
monitor repro_Transect_precision, thin(5)
monitor reproyesno_Transect_precision, thin(5)
monitor surv_Transect_precision, thin(5)
monitor grwth_Transect_precision, thin(5)
update 50000
parameters to "out2.Rdump", chain(1)
coda *, stem(sim.2/CODA)
samplers to sim.2/samplers.csv
update 0
model clear
exit
