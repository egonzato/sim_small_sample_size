using CSV, LSurvival, GLM, Plots, StatsModels, DataFrames, Distributions, Random
# predict
function predictsurv(model,data)
    ch=cumsum(model.bh[:,1])
    β=coef(model)
    X=Matrix(data[:,coefnames(model)])
    η=X * β
    st= exp.( .- ch' .* exp.(η))
    return st
end
# ipw
function ipw(data::DataFrame, exposure::Symbol, confounders::Vector{Symbol}, type::String)
    # Construct the formula string
    conf = join(string.(confounders), " + ")
    fmla = @eval @formula($(exposure) ~ $(Meta.parse(conf)))
    # Fit the logistic regression model
    model = glm(fmla, data, Binomial(), LogitLink())
    # Get predicted probabilities
    prob = predict(model,data,type=:response)
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
# gcomp
function gcomp(
    data::DataFrame,
    time::Symbol,
    status::Symbol,
    treatment::Symbol,
    confounders::Vector{Symbol};
    bootstrap::Bool = false,
    B::Int = 1000,
    rng::AbstractRNG = Random.GLOBAL_RNG
)
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
