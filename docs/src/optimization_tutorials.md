# Non-Allocating Forward-Mode L2 Collocation Loss

The following is an example of a loss function over the collocation that
is non-allocating and compatible with forward-mode automatic differentiation:

```julia
using PreallocationTools
du = PreallocationTools.dualcache(similar(prob.u0))
preview_est_sol = [@view estimated_solution[:, i] for i in 1:size(estimated_solution, 2)]
preview_est_deriv = [@view estimated_derivative[:, i]
                     for i in 1:size(estimated_solution, 2)]

function construct_iip_cost_function(f, du, preview_est_sol, preview_est_deriv, tpoints)
    function (p)
        _du = PreallocationTools.get_tmp(du, p)
        vecdu = vec(_du)
        cost = zero(first(p))
        for i in 1:length(preview_est_sol)
            est_sol = preview_est_sol[i]
            f(_du, est_sol, p, tpoints[i])
            vecdu .= vec(preview_est_deriv[i]) .- vec(_du)
            cost += sum(abs2, vecdu)
        end
        sqrt(cost)
    end
end
cost_function = construct_iip_cost_function(
    f, du, preview_est_sol, preview_est_deriv, tpoints)
```