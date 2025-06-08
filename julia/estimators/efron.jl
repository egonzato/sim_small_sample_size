using DataFrames, StatsBase, Plots

# Assume df is your DataFrame with columns: :time, :y (1=event, 0=censored), :trt (or other covariates)
# Assume you have fit your Cox model and have coefficients:
beta = coef(model)  # e.g., [β_trt]

# Compute linear predictor and risk score for each subject
df.eta = df.trt .* beta[1]  # Extend if you have more covariates
df.risk_score = exp.(df.eta)

# Get unique event times (only where events occurred)
event_times = sort(unique(df.time[df.y .== 1]))
n_events = length(event_times)

# Prepare storage
baseline_hazard = zeros(n_events)
cumulative_hazard = zeros(n_events)

for (j, t) in enumerate(event_times)
    # Risk set: subjects at risk just before time t
    risk_set = (df.time .>= t)
    # Events at time t
    event_at_t = (df.time .== t) .& (df.y .== 1)
    d_j = sum(event_at_t)  # Number of events at t
    # Risk scores
    sum_risk = sum(df.risk_score[risk_set])
    sum_event_risk = sum(df.risk_score[event_at_t])
    # Efron estimator for ties
    h0 = 0.0
    for l in 0:(d_j-1)
        h0 += 1.0 / (sum_risk - l/d_j * sum_event_risk)
    end
    baseline_hazard[j] = h0
    cumulative_hazard[j] = sum(baseline_hazard[1:j])
end

# Calculate survival probabilities for a covariate pattern, e.g., trt=1
eta_pred = 1 * beta[1]  # adjust for your covariate(s)
surv_prob = exp.(-cumulative_hazard .* exp(eta_pred))

# Create DataFrame with results
results = DataFrame(
    time = event_times,
    baseline_hazard = baseline_hazard,
    cumulative_hazard = cumulative_hazard,
    survival_prob = surv_prob
)
# load comparison dataset
compare=CSV.read("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\compare.csv", DataFrame)
# Plot the survival curve
plot(event_times, surv_prob,
     xlabel="Time",
     ylabel="Predicted Survival Probability",
     title="Manual Efron Survival Curve",
     legend=false,
     lw=2)
plot!(compare.time,compare.surv,
lw=1.5)
