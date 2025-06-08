# loading libraries
using Random, Distributions, Plots, Statistics, StatsModels, DataFrames, StatsPlots, CSV
# sample size
n=1000
# create vector for age
age=rand(Normal(45,3),n)
eps=rand(Uniform(-0.12,0.12),n)
# probability of treatment
logittrt=0.015.*age.+eps
ptrt=(1 ./(1 .+exp.(-logittrt)))
trt=rand.(Bernoulli.(ptrt))
histogram(ptrt)
# same for the outcome
logity = 0.05 .+ 0.007 .* age .- 0.25 .* trt .+ eps
py=(1 ./(1 .+exp.(-logity)))
histogram(py)
y=rand.(Bernoulli.(py))
# simulate time
time=-log.(rand(Uniform(0,1),n)) ./ (exp.(.- 0.2 .* trt .+ 0.002 .* age .+ rand(Normal(0.5,1),n)))
histogram(time)
# create df
df=DataFrame(age=age, trt=trt, time=time, y=y)
df.weight = rand(Uniform(0.5, 1.5), nrow(df))  # Example weights, replace with actual weights if available
# export df
CSV.write("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", df)