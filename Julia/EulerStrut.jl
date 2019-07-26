using BoundaryValueDiffEq, OrdinaryDiffEq, Plots

# Initialise parameters
P = 2*pi + 0.1 # Kickstart the problem
u0 = [0.0, 1.0*P] # Initial guess
tspan = (0.0, 1.0) # Domain

# Define ODE function
function Odes!(du, u, p, t)
    theta = u[1]
    phi = u[2]
    du[1] = phi
    du[2] = -P^2*sin(theta)
end 

# Define BC function
function Bc!(residual, u, p, t)
    residual[1] = u(0)[1]
    residual[2] = u(1)[1]
end

# First solve it in terms of shooting 
bvp = BVProblem(Odes!, Bc!,  u0, tspan)
@time sol = solve(bvp, Shooting(Vern7()))

plot(sol, vars=(1))

# # Now try it using collocation
# bvp2 = TwoPointBVProblem(Odes!, Bc!, [sin(P*t), P*cos(P*t)], tspan)
# sol2 = solve(bvp2, GeneralMIRK4(), dt=0.05)

# plot(sol2)