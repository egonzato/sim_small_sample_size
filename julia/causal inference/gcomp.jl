using CSV, DataFrames, LSurvival, StatsModels, StatsBase, CategoricalArrays, Distributions, Plots, Bootstrap
# Load dataset
df=CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv", DataFrame)
# Add enter column
df.enter =zeros(length(df.y))
# Ensure correct types
df.trt = Int.(df.trt)   # Or categorical(df.trt) if it's a group
df.y = Int.(df.y)
# Fit Cox model
model = LSurvival.coxph(@formula(Surv(enter, time, y) ~ trt + age), df,
#                                 wts=df.weight,
                                 ties = "efron")
# treated and untreated df
ctrl = copy(df)
trt = copy(df)
ctrl.trt .= 0  # Use .= for in-place modification
trt.trt .= 1
# predict
pred_ctrl=predictsurv(model,ctrl); pred_trt=predictsurv(model,trt)
# store in dataset
df.pred_ctrl=pred_ctrl[:,size(pred_ctrl)[2]]; df.pred_trt=pred_trt[:,size(pred_trt)[2]]
# now calculate ate
ate=log.((.- log.(mean(df.pred_trt))) ./ (.- log.(mean(df.pred_ctrl))))
# creat nboot and ates vector
nobs=size(df)[1]
ates=Float64[]
# bootstrap
for i in 1:1000
    newdf=df[sample(1:nobs,nobs;replace=true),:]
    model = LSurvival.coxph(@formula(Surv(enter, time, y) ~ trt + age), newdf)
    # treated and untreated df
    ctrl = copy(newdf)
    trt = copy(newdf)
    ctrl.trt .= 0  # Use .= for in-place modification
    trt.trt .= 1
    # predict
    pred_ctrl=predictsurv(model,ctrl); pred_trt=predictsurv(model,trt)
    # store in dataset
    newdf.pred_ctrl=pred_ctrl[:,end]; newdf.pred_trt=pred_trt[:,end]
    # now calculate ate
    ate=log.((.- log.(mean(newdf.pred_trt))) ./ (.- log.(mean(newdf.pred_ctrl))))
    push!(ates,ate)
end
# define function
using DataFrames, StatsBase, LSurvival, StatsModels
# function
function gcompstd(dts, i)
    # shuffled dataset
    data = dts[i, :]
    # create enter indicator
    data[!, :enter] = zeros(nrow(data))  
    # modify boolean to dummy
    data[!, :trt] = Int.(data[!, :trt])
    data[!, :y] = Int.(data[!, :y])
    # fit model
    model = LSurvival.coxph(@formula(Surv(enter, time, y) ~ trt + age), data)
    # create counterfactucal dataset
    ctrl=trt=deepcopy(data)
    ctrl[!, :trt] .= 0
    trt[!, :trt] .= 1
    # predict
    pred_ctrl, pred_trt=predict(model, ctrl), predict(model, trt)
    # store predictions
    data[!, :pred_ctrl] = pred_ctrl
    data[!, :pred_trt] = pred_trt
    # calculate ate
    ate = log(-log(mean(pred_trt)) / -log(mean(pred_ctrl)))
    return ate
end
# Example bootstrap
ates_boot = Bootstrap.bootstrap(gcompstd, df, BasicSampling(1000))
function gcomp(
    data::DataFrame,
    time::Symbol,
    status::Symbol,
    treatment::Symbol,
    confounders::Vector{Symbol};
    bootstrap::Bool = false,
    B::Int = 1000,
    rng::AbstractRNG = Random.GLOBAL_RNG)
    # Build formula
    rhs_vars = [treatment; confounders...]
    rhs = join(string.(rhs_vars), " + ")
    formula_str = "@formula(Surv($(string(time)), $(string(status))) ~ $rhs)"
    f = eval(Meta.parse(formula_str))
    model = coxph(f, data)
    # Predict for all treated and all controls
    ctrl = deepcopy(data)
    trt = deepcopy(data)
    ctrl[!, treatment] .= 0
    trt[!, treatment] .= 1
    pred_ctrl = predictsurv(model, ctrl)[:, end]
    pred_trt = predictsurv(model, trt)[:, end]
    ate = log.((.- log.(mean(pred_trt))) ./ (.- log.(mean(pred_ctrl))))

    # Bootstrap if requested
    if bootstrap
        n = nrow(data)
        ates = Vector{Float64}(undef, B)
        for b in 1:B
            idxs = rand(rng, 1:n, n)  # sample with replacement
            data_b = data[idxs, :]
            # Fit model and compute ATE on bootstrap sample
            model_b = coxph(f, data_b)
            ctrl_b = deepcopy(data_b)
            trt_b = deepcopy(data_b)
            ctrl_b[!, treatment] .= 0
            trt_b[!, treatment] .= 1
            pred_ctrl_b = predictsurv(model_b, ctrl_b)[:, end]
            pred_trt_b = predictsurv(model_b, trt_b)[:, end]
            ates[b] = log.((.- log.(mean(pred_trt_b))) ./ (.- log.(mean(pred_ctrl_b))))
        end
        std_boot = std(ates)
        return (ate=ate, std=std_boot)
    else
        return (ate=ate)
    end
end
