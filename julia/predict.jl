# libraries
using CSV, LSurvival, DataFrames, Distributions
# write function
function predict(model,data)
    ch=cumsum(model.bh[:,1])
    β=coef(model)
    X=Matrix(data[:,coefnames(model)])
    η=X * β
    st= exp.( .- ch' .* exp.(η))
    return st
end
# read dataset
df = CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", DataFrame)
# compute model
# Add enter column
df.enter =zeros(length(df.y))
# Ensure correct types
df.trt = Int.(df.trt)   # Or categorical(df.trt) if it's a group
df.y = Int.(df.y)
df.weight = rand(Uniform(0.5, 1.5), nrow(df))  # Example weights, replace with actual weights if available
# Fit Cox model
model = coxph(@formula(Surv(enter, time, y) ~ trt + age), df,ties = "efron")
# use fcuntion
st1=predict(model,df)
# last row for ate 
st1ate=st1[:,size(st1)[2]]