using Random, Distributions, LSurvival, StatsModels, StatsPlots
# number of observations
n=25000
# simulate big population from which the loop will sample
age=rand(Normal(55,15),n)
bmi=rand(Normal(25,5),n)
sex=rand(Bernoulli(0.5),n)
# simulate treatment
logittrt = -11 .+ 0.05.*sex .+ 0.19.*bmi .+ 0.1.*age .+ rand(Uniform(-0.8,0.8), n)
ptrt=(1 ./(1 .+exp.(-logittrt)))
trt=rand.(Bernoulli.(ptrt))
histogram(ptrt)
# simulate event
logitevent = -2.7 .+ 0.04 .* sex .+ 0.07 .* bmi .+ 0.012 .* age .- 0.6 .* trt .+ rand(Uniform(-0.7,0.7), n)
pevent = 1 ./ (1 .+ exp.(-logitevent))
event=rand.(Bernoulli.(pevent))
histogram(pevent)
# simulate survival time
time = -log.(rand(Uniform(0, 1), n)) ./ exp.(-0.2 .* trt .+ 0.002 .* age .+ 0.03 .* bmi .+ 0.05 .* sex .+ rand(Normal(0.5, 1), n))
# create df
population=DataFrame(time=time,event=event,age=age,bmi=bmi,sex=sex,trt=trt)
# distribution of age in two different treatment groups
histogram(bmi, group=trt, legend=:topleft, alpha=0.6, bins=30)
# write to a dataframe