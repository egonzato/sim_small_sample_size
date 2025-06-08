using CSV, LSurvival, GLM, Plots, StatsModels, DataFrames
# load dataset
df = CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", DataFrame)
# logistic regression of treatment against age
trt_model=glm(@formula(trt~age),df,Binomial(),LogitLink())
# now predict on the same data
df.prob=predict(trt_model,rwd,type=:response)
# type of weights
type="stabilized"
# calculate weights
df.weight=ifelse.(type .=="normal",
                 (df.trt ./ df.prob) .+((1 .- df.trt) ./ (1 .- df.prob)),
                 ifelse.(df.trt .==0,
                        (1 .- mean(df.trt)) ./ (1 .- df.prob),
                        mean(df.trt) ./ df.prob))
# distribution of weights in treated and untreated
histogram(df.weight[df.trt .== 1])
histogram!(df.weight[df.trt .== 0])
# now the cox model
df.enter =zeros(length(df.y))
# Ensure correct types
df.trt = Int.(df.trt)   # Or categorical(df.trt) if it's a group
df.y = Int.(df.y)
# Fit Cox model
model = LSurvival.coxph(@formula(Surv(enter, time, y) ~ trt),
                        df,
                        wts=df.weight,
                        ties = "efron")
# function
function ipw(data::DataFrame, exposure::Symbol, confounders::Vector{Symbol}, type::String)
    # Construct the formula string
    conf = join(string.(confounders), " + ")
    fmla = @eval @formula($(exposure) ~ $(Meta.parse(conf)))
    # Fit the logistic regression model
    model = glm(fmla, data, Binomial(), LogitLink())
    # Get predicted probabilities
    prob = predict(model)
    # Marginal probability of treatment
    p = mean(data[!, exposure])
    # Calculate weights
    trt = data[!, exposure]
    weights = ifelse.(type == "normal",
        (trt ./ prob) .+ ((1 .- trt) ./ (1 .- prob)),
        ifelse.(trt .== 1,
            p ./ prob,
            (1 - p) ./ (1 .- prob)
        )
    )
    return weights
end
# try function on dataset
weights = ipw(df, :trt, [:age], "stabilized")

