# download MatchIt form GitHub
# Pkg.add(url="http://github.com/eohne/MatchIt.jl")
# ?MatchIt
# import libraries
using CSV, DataFrames, LSurvival, StatsModels, Plots, MatchIt
# Load dataset
df = CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", DataFrame)
# perform propensity score matching 2:1 ratio
