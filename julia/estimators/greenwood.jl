using LinearAlgebra

# Assume you have: event_times, baseline_hazard, cumulative_hazard, df.risk_score, df.eta, etc.

# Greenwood variance for cumulative hazard at each event time
greenwood_var = zeros(length(event_times))
for j in 1:length(event_times)
    v = 0.0
    for k in 1:j
        t_k = event_times[k]
        risk_set = df.time .>= t_k
        d_k = sum((df.time .== t_k) .& (df.y .== 1))
        sum_risk = sum(df.risk_score[risk_set])
        v += d_k / (sum_risk^2)
    end
    greenwood_var[j] = v
end

# For a given eta_pred (e.g., trt=1)
eta_pred = 1 * beta[1]
surv_prob = exp.(-cumulative_hazard .* exp(eta_pred))

# Variance of survival
var_surv = (surv_prob.^2) .* greenwood_var .* exp(2 * eta_pred)

# 95% confidence intervals (log-log scale, as in R)
loglog_surv = log.(-log.(surv_prob))
se_loglog = sqrt.(greenwood_var) * exp(eta_pred)
z = 1.96  # for 95% CI

upper = exp.(-exp.(loglog_surv - z * se_loglog))
lower = exp.(-exp.(loglog_surv + z * se_loglog))

# Plot
plot(event_times, surv_prob, label="Survival", lw=2)
plot!(event_times, upper, label="Upper 95% CI", ls=:dash)
plot!(event_times, lower, label="Lower 95% CI", ls=:dash)
xlabel!("Time")
ylabel!("Predicted Survival Probability")
title!("Predicted Survival Curve with 95% CI")
