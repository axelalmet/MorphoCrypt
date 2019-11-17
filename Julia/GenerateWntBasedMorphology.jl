using BoundaryValueDiffEq, OrdinaryDiffEq
using Plots

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

# Define the Wnt function
W = S-> exp.(-(S./sigma).^2)


# Define the initial growth
sOld = S-> S
gammaOld = S-> 1.0

# Define the initial foundation
px = S-> S
py = S-> 0.0

# Initial guess obtained from AUTO-07p for linear elastic foundation
initValues = [0.0, 0.0, -0.3351, -98.07, 0.0, 0.0, 0.9871]
SSpan = (0.0, 0.5) # Solution span

function RodOdes!(du, u, p, S)
    # Define solution components
    s, x, y, nx, ny, theta, m = u

    gammaNew(S) = gammaOld(S).*(1 + 0.5*dt.*W(sOld(S)))./(1 - 0.5*dt*W(s))

    du[1] = gammaNew(S) # s' = gamma
    du[2] = gammaNew(S).*cos(theta) # x' = gamma * cos(theta)
    du[3] = gammaNew(S).*sin(theta) # y' = gamma * sin(theta)
    du[4] = K.*gammaNew(S).*(x - px(S)) # nx' = k*gamma*(x - px)
    du[5] = K.*gammaNew(S).*(y - py(S)) # ny' = k*gamma*(y - py)
    du[6] = gammaNew(S).*m # theta' = gamma*m
    du[7] = gammaNew(S).*(nx.*sin(theta) - ny.*cos(theta)) # m' = gamma(nx*sin(theta) - ny*cos(theta))

end

function Bcs!(residual, u, p, S)

    residual[1] = u(0.0)[1] # s(0) = 0
    residual[2] = u(0.0)[2] # x(0) = 0
    residual[3] = u(0.0)[5] # ny(0) = 0
    residual[4] = u(0.0)[6] # theta(0) = 0
    residual[5] = u(0.5)[2] - 0.5 # x(L) = L/2
    residual[6] = u(0.5)[3] # y(L) = 0
    residual[7] = u(0.5)[6] # theta(L) = 0

end

# # Get the initial solution. We use the shooting method.
bvp = BVProblem(RodOdes!, Bcs!, initValues, SSpan)
@time sol = solve(bvp, Shooting(Vern7(lazy=false)))

# plot(sol, vars=(2, 3))

# Update growth
sNew(S) = sol(S, idxs=1)
# print(gammaOld(0)*(1 + 0.5*dt.*W(sNew(0)))./(1 - 0.5*dt.*W(sOld(0))))
gammaNew(S) = gammaOld(S).*(1 + 0.5*dt.*W(sNew(S)))./(1 - 0.5*dt.*W(sOld(S)))

print(gammaNew([0.0]))

# gammaOld(S) = gammaNew(S)

# # Update the foundation. We need the solutions for x and y to do this.
# xNew = S -> sol(S, idxs=2)
# yNew = S -> sol(S, idxs=3) 
# pxOld = px(S)
# pyOld = S -> py(S)

# px = S -> (pxOld(S) + 0.5*dt*rho*(xNew(S) - xOld(S) + pxOld(S)))./(1 + 0.5*dt*rho)
# py = S -> (pyOld(S) + 0.5*dt*rho*(yNew(S) - yOld(S) + pyOld(S)))./(1 + 0.5*dt*rho)

# print(px[0 0.1 0.2])
# print(py[0 0.1 0.2])

# # Update the solutions and current arc length
# solOld = sol
# sOld = S -> solOld(S, idxs=1)

#  # Now try to solve it again
# initValues = solOld(0)
# bvpNew = BVProblem(RodOdes!, Bcs!, initValues, SSpan)
# @time solNew = solve(bvp, Shooting(Vern7()))