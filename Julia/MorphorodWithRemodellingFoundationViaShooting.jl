using BoundaryValueDiffEq, OrdinaryDiffEq, Plots

# Define the parameters
dt = 0.05 # time step 
g = 1.0 # Growth rate
kf = 0.01 # ratio of rod stiffness to foundation stiffness
h = 15.0
w = 10.0
L0 = 125.0
sigma = 0.16 # Width of wnt gradient

# Define the foundation stiffness
K = 12.0*kf*L0^4/(w*h^3)
rho = 10.0 # relaxation rate of foundation

# Rod stiffness (homogeneous this time)
Eb = 1.0

# Define the Wnt function
W = S -> exp(-(S./sigma).^2)

# Define the initial growth
gamma = 1 + dt

# Define the initial foundation
local px = S -> S
local py = S -> 0.0

# Initial guess for linear elastic foundation
initValues = [0.0, 0.0, -0.3351, -98.07, 0.0, 0.0, 0.9871]
SSpan = (0.0, 0.5) # Solution span

function RodOdes!(du, u, p, S)
    # Define solution components
    s, x, y, nx, ny, theta, m = u

    # gamma = S -> gammaOld(S).*(1 + 0.5*dt*W(sOld))./(1 - 0.5*dt*W(s))
    du[1] = gamma # s' = gamma
    du[2] = gamma.*cos(theta) # x' = gamma * cos(theta)
    du[3] = gamma.*sin(theta) # y' = gamma * sin(theta)
    du[4] = K.*gamma.*(x - px(S)) # nx' = k*gamma*(x - px)
    du[5] = K.*gamma.*(y - py(S)) # ny' = k*gamma*(y - py)
    du[6] = gamma.*m./Eb # theta' = gamma/Eb*m
    du[7] = gamma.*(nx.*sin(theta) - ny.*cos(theta)) # m' = gamma(nx*sin(theta) - ny*cos(theta))
end

function Bcs!(residual, u, p, S)

    residual[1] = u(0)[1] # s(0) = 0
    residual[2] = u(0)[2] # x(0) = 0
    residual[3] = u(0)[5] # ny(0) = 0
    residual[4] = u(0)[6] # theta(0) = 0
    residual[5] = u(0.5)[2] - 0.5 # x(L) = L/2
    residual[6] = u(0.5)[3] # y(L) = 0
    residual[7] = u(0.5)[6] # theta(L) = 0

end

# # Get the initial solution. We use the shooting method.
bvp = BVProblem(RodOdes!, Bcs!, initValues, SSpan)
@time sol = solve(bvp, Shooting(Vern7()))

# plot(sol, vars=(2, 3))

# print(sol(0, idx=1:7))    
# # Update the foundation. We need the solutions for x and y to do this.
xOld= S -> sol(S, idxs=2)
yOld = S -> sol(S, idxs=3) 

 # S->  (pxOld(S) + 0.5*dt*rho*(xNew(S) - xOld(S) - pxOld(S)))./(1 + 0.5*dt*rho)

local pxOld = px
local py = S->  (pxOld(S) + dt*rho*(xOld(S) - pxOld(S)))

# S-> (pyOld(S) + 0.5*dt*rho*(yNew(S) - yOld(S) - pyOld(S)))./(1 + 0.5*dt*rho)
local pyOld = py
local py = S-> (pyOld(S) + dt*rho*(yOld(S) - pyOld(S)))

print(px([0.1, 0.2]))
print(py([0.1, 0.2]))

# # Update the solutions and current arc length
solOld = sol
# sOld = S -> solOld(S, idx=1)

# # Now try to solve it again
# initValues = solOld(0)
# bvpNew = BVProblem(RodOdes!, Bcs!, initValues, SSpan)
# solNew = solve(bvp, Shooting(Vern7()))