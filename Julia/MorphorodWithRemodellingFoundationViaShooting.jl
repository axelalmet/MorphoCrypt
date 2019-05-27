using DifferentialEquations
using Interpolations

# Import initial guess
initSolData = readdlm("../Data/planarmorphorodsk0p02L29_sol_1")
solMesh = initSolData[:,1]'
initS = initSolData[:,2]'
initX = initSolData[:,3]'
initY = initSolData[:,4]'
initF = initSolData[:,5]'
initG = initSolData[:,6]'
initTheta = initSolData[:,7]'
initM = initSolData[:,8]'

function Odes!(du, u, p, t)

    X = u[1,:]
    Y = u[2,:]
    Nx = u[3,:]
    Ny = u[4,:]
    theta = u[5,:]
    m = u[6,:]

    dx = gamma.*cos(Theta)
    dy = gamma.*sin(Theta)
    dnx = K.*gamma.*(X - px)
    dny = K.*gamma.*(Y - py)
    dtheta = gamma.*m.*Eb.^(-1)
    dm = gamma.*(Nx.*sin(Theta) - Ny.*cos(Theta))

    du = [dx; dy; dnx; dny; dtheta; dm]
end

function Bcs!(residual, u, p, t)

    residual[1] = u[end][2] - L # x(L) = L
    residual[2] = u[end][3] # y(L) = 0
    residual[3] = u[end][5] # theta(L) = 0

end