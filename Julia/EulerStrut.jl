using BoundaryValueDiffEq, OrdinaryDiffEq, Plots

# Initialise parameters
P = 2*pi + 0.1 # Kickstart the problem
u0 = [0.0, 1.0*P] # Initial guess
tspan = (0.0, 1.0) # Domain

# Define ODE function. We write it in this way to try and be more efficient with the code
function Odes!(du, u, p, t)
    du[1] = u[2] # theta' = phi
    du[2] = -P^2*sin(u[1]) # phi' = -P^2*sin(theta) 
end 

# Define BC function
function Bc!(residual, u, p, t)
    residual[1] = u(0)[1]
    residual[2] = u(1)[1]
end

# First solve it in terms of shooting 
bvp = BVProblem(Odes!, Bc!,  u0, tspan)
@time sol = solve(bvp, Vern7())

sol(0)