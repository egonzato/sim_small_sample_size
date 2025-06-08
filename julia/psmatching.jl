# libraries
using CSV, StatsModels, GLM, Random, Distributions, DataFrames, StatsPlots
weights=rand(Uniform(0,1),1000)
trt=rand.(Bernoulli.(weights))
# create weights, treatment and bind into dataset
df=DataFrame(weights=weights,trt=trt)
# filter
untrt=df[df.trt .==false,:]
trt=df[df.trt .==true,:]
# set indicator
i=1
# calculate difference
diff=abs.(trt.weights[i] .- untrt.weights)
# find index of minimum
val, idx=findmin(diff)
# create new df with the new observations and the cluster identificator
clustered = DataFrame(weights=Float64[], cluster=Int[])
push!(clustered, (trt.weights[i],idx))
push!(clustered, (untrt.weights[idx],idx))
delete!(untrt,idx)
# try on all observations
df=DataFrame(weights=weights,trt=trt)
# filter df
untrt=df[df.trt .==false,:]
trt=df[df.trt .==true,:]
# empty dataset where to store clustered obs
clustered = DataFrame(weights=Float64[], cluster=Int[])
# loop through
for i in 1:size(trt)[1]
    if nrow(untrt) == 0
        break  
    end
    diff=abs.(trt.weights[i] .- untrt.weights)
    # find index of minimum
    val, idx=findmin(diff)
    # create new df with the new observations and the cluster identificator
    push!(clustered, (trt.weights[i],idx))
    push!(clustered, (untrt.weights[idx],idx))
    delete!(untrt,idx)
end
# try on dataset
rwd=CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", DataFrame)
# distribution of age under treated untreated
@df rwd histogram(:age, group = :trt, legend = :topleft, alpha = 0.6, bins=30)
# model
model=glm(@formula(trt~age),rwd,Binomial(),LogitLink())
# predict weights
rwd.wi=predict(model,rwd,type=:response)
# create empy df
untrt=rwd[rwd.trt .==false,:]
trt=rwd[rwd.trt .==true,:]
# empty dataset where to store clustered obs
clustered = DataFrame()
# loop through
maxiter=min(nrow(trt),nrow(untrt))
for i in 1:maxiter
    if nrow(untrt) == 0
        break  
    end
    diff=abs.(trt.wi[i] .- untrt.wi)
    # find index of minimum
    val, idx=findmin(diff)
    # create new df with the new observations and the cluster identificator
    push!(clustered, merge(NamedTuple(trt[i, :]), (cluster=i,)))
    push!(clustered, merge(NamedTuple(untrt[idx, :]), (cluster=i,)))
    delete!(untrt,idx)
end
# distribution of age in the matched dataset
@df clustered histogram(:age, group = :trt, legend = :topleft, alpha = 0.6, bins=30)
