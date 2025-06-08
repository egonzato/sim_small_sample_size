using Survival, DataFrames, CSV, Plots, StatsModels
# load dataset
df = CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in sample sizes\\Data\\julia\\df.csv", DataFrame)
# create event variable
df.event=EventTime.(df.time, df.y .== 1)
# fit km
km = fit(KaplanMeier, df.time, df.y)
ci = confint(km)
lower = [x[1] for x in ci]
upper = [x[2] for x in ci]
# Create basic plot
plot(km.events.time, km.survival, 
    label="KM Curve", 
    xlabel="Time", 
    ylabel="Survival Probability",
    linewidth=2)
plot!(km.events.time, lower, fillrange=upper, 
    fillalpha=0.2, 
    label="95% CI",
    legend=:bottomleft)
# plot the two treatment
# Split data by treatment
df_trt1 = filter(row -> row.trt == true, df)
df_trt0 = filter(row -> row.trt == false, df)
# Fit separate KM estimators
km_trt1 = fit(KaplanMeier, df_trt1.time, df_trt1.y)
km_trt0 = fit(KaplanMeier, df_trt0.time, df_trt0.y)
# Create comparative plot
plot(km_trt1.events.time, km_trt1.survival, 
    label="Treatment Group", 
    xlabel="Time", 
    ylabel="Survival Probability",
    linewidth=2)
plot!(km_trt0.events.time, km_trt0.survival, 
    label="Control Group",
    linewidth=2)
# cox model
model=Survival.coxph(@formula(event~trt), df)
# plot(survivalfit((@formula(Surv(time, y) ~ 1)), df))
coef(model)
vcov(model)
# content of the model