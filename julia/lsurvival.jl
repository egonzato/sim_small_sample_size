using CSV, DataFrames, LSurvival, StatsModels, StatsBase, CategoricalArrays, Distributions, Plots
# Load dataset
df = CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", DataFrame)
# Add enter column
df.enter =zeros(length(df.y))
# Ensure correct types
df.trt = Int.(df.trt)   # Or categorical(df.trt) if it's a group
df.y = Int.(df.y)
df.weight = rand(Uniform(0.5, 1.5), nrow(df))  # Example weights, replace with actual weights if available
# Fit Cox model
model = LSurvival.coxph(@formula(Surv(enter, time, y) ~ trt), df,
#wts=df.weight,
ties = "efron")
# baseline hazard
ch=cumsum(model.bh[:,1])
bh=model.bh[:,1]
tevents=model.bh[:,4]
coeff=repeat(coef(model),size(model.bh)[1])
trt=fill(1,size(model.bh)[1])
# baseline from cumulative
preds=DataFrame(ch=ch,bh=bh,tevents=tevents,coeff=coeff,trt=trt)
# calculate probabilities at each time
preds.HR=exp.(preds.trt .* preds.coeff)
preds.preds=(exp.(-preds.ch)).^preds.HR
# now plot
plot(preds.tevents, preds.preds,
     xlabel="Time",
     ylabel="Predicted Survival Probability",
     title="Predicted Survival Curve",
     legend=false,
     lw=2)