using BoundaryValueDiffEq, OrdinaryDiffEq, Plots

# Initialise parameters
P = 4.75 # Kickstart the problem
L = pi
u0 = [0.0, 0.0, -6.5184, 0.0, 0.0, 3.5899] # Initial guess
tspan = (0.0, 1.0) # Domain

# Define ODE function before self-contact
function Odes!(du, u, p, t)
    x = u[1]
    y = u[2]
    Q = u[3]
    N = u[4]
    theta = u[5]
    k = u[6]

    du[1] = L*cos(theta)
    du[2] = L*sin(theta)
    du[3] = -L*k*N
    du[4] = L*(k*Q + P)
    du[5] = L*k
    du[6] = L*N

end 

# Define ODE function before self-contact
function SelfContactOdes!(du, u, p, t)
    x = u[1]
    y = u[2]
    Q = u[3]
    N = u[4]
    theta = u[5]
    k = u[6]

    du[1] = L*cos(theta)
    du[2] = L*sin(theta)
    du[3] = -L*k*N
    du[4] = L*(k*Q + P)
    du[5] = L*k
    du[6] = L*N

end 

# Define BC function before self-contact
function PreContactBc!(residual, u, p, t)
    residual[1] = u(0)[1]
    residual[2] = u(0)[2]
    residual[3] = u(0)[4]
    residual[4] = u(1)[4]
    residual[5] = u(0)[5]
    residual[6] = u(1)[5] - L
end

# Define the BCs during self-contact at a point
function SelfPointContactBc!(residual, u, p, t)
    residual[1] = u(0)[1]
    residual[2] = u(0)[2]
    residual[3] = u(0.5)[1]
    residual[4] = u(0.5)[5] - 0.5*pi
    residual[5] = u(1)[1]
    residual[6] = u(1)[5] - L
end

# Define the BCs during self-contact along a region
function SelfBeforeRegionContactBc!(residual, u, p, t)
    residual[1] = u(0)[1]
    residual[2] = u(0)[2]
    residual[3] = u(0)[5]
    residual[4] = u(0.5)[1]
    residual[5] = u(0.5)[5] - 0.5*pi
    residual[6] = u(0.5)[6]
end

# Define the BCs during self-contact along a region
function SelfAfterRegionContactBc!(residual, u, p, t)
    residual[1] = u(0)[1]
    residual[2] = u(0)[5] - 0.5*pi
    residual[3] = u(0)[6]
    residual[4] = u(0)[3] - Qsc
    residual[5] = u(0.5)[1]
    residual[6] = u(0.5)[5] - L
end

# First solve it in terms of shooting 
precontactbvp = BVProblem(Odes!, PreContactBc!,  u0, tspan)
@time sol = solve(precontactbvp, Shooting(Vern7()), force_dtmin=true)

p1 = plot(sol, vars=(1, 2), reuse=false)

# Let's try solving the self-contact solution
P = 5.247

u0 = [0, 0, -7.1859, 0, 0, 3.9207]

selfpointcontactbvp = BVProblem(Odes!, SelfPointContactBc!,  u0, tspan)
@time sol = solve(selfpointcontactbvp, Shooting(Vern7()), force_dtmin=true)

p2 = plot(sol, vars=(1, 2), reuse=false)

# Self-contact along a region
P = 10.4
u0 = [0, 0, -10.5603, 0, 0, 4.5957]
tspan = (0, 0.5)

selfbeforeregioncontactbvp = BVProblem(Odes!, SelfBeforeRegionContactBc!,  u0, tspan)
@time sol1 = solve(selfbeforeregioncontactbvp, Shooting(Vern7()), force_dtmin=true)

u0 = sol1[:, end] + [0, 0.1, 0, 0, 0, 0]
Qsc = sol1[end,3]
selfafterregioncontactbvp = BVProblem(Odes!, SelfAfterRegionContactBc!,  u0, tspan)
@time sol2 = solve(selfafterregioncontactbvp, Shooting(Vern7()), force_dtmin=true)

p3 = plot(sol1, vars=(1, 2))
p3 = plot!(sol2, vars=(1, 2))

plot(p1, p2, p3)
