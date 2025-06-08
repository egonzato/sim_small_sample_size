using CSV, DataFrames, LSurvival, StatsModels, StatsBase, CategoricalArrays, LinearAlgebra, Distributions, Plots
# Load dataset
df = CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", DataFrame)
# Add enter column
df.enter =zeros(length(df.y))
# Ensure correct types
df.trt = Int.(df.trt)   # Or categorical(df.trt) if it's a group
df.y = Int.(df.y)
# Fit Cox model
model = LSurvival.coxph(@formula(Surv(enter, time, y) ~ trt + age), df,
wts=df.weight,
ties = "efron")
# extract coefficients
β=coef(model)
# cumulative hazard
ch=cumsum(model.bh[:,1])
tevents=model.bh[:,4]
# new matrix
X=Matrix(df[:,[coefnames(model)]])
# linear predictor
η=X * β
# survival probability
st= exp.( .- ch' .* exp.(η))
# extract last row
atest=st[:,size(st)[2]]
# now plot
plot(tevents,st[253, :],
     xlabel="Time",
     ylabel="Predicted Survival Probability",
     title="Predicted Survival Curve",
     legend=false,
     lw=2)

# Plot all survival curves
plt = plot()  # Start with an empty plot
for i in 1:size(st)[1]
    plot!(plt, tevents, st[i, :], color=:auto, alpha=0.5, label=false)
end

xlabel!(plt, "Time")
ylabel!(plt, "Predicted Survival Probability")
title!(plt, "Individual Survival Curves")
display(plt)