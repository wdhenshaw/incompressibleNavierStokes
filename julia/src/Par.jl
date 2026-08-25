
# BC type codes
const PERIODIC        = -1
const DIRICHLET       =  1
const NO_SLIP_WALL    =  2
const SLIP_WALL       =  3
const INFLOW          =  4
const OUTFLOW         =  5
const PRESSURE_INFLOW =  6

# component codes
const UC = 1
const VC = 2
const PC = 3

mutable struct Par
    # --- physics ---
    nd::Int
    nu::Float64
    rho::Float64
    # --- domain ---
    xa::Float64; xb::Float64
    ya::Float64; yb::Float64
    # --- grid ---
    N0::Int
    Nx::Int; Ny::Int
    Ngx::Int; Ngy::Int
    dx::Float64; dy::Float64
    dr::Vector{Float64}          # dr[1]=dx in r, dr[2]=dy in s
    numGhost::Int
    orderInSpace::Int
    gid::Matrix{Int}             # gid[side,axis]: boundary indices  (2×2)
    dim::Matrix{Int}             # dim[side,axis]: total array bounds (2×2)
    isCartesian::Bool
    x::Array{Float64,3}          # x[i1,i2,1:2] physical coords
    rx::Array{Float64,4}         # rx[i1,i2,m1,m2] metric: dr_m/dx_n (curvilinear)
    map::String
    # --- time stepping ---
    tf::Float64
    tp::Float64
    tFinal::Float64
    dtMax::Float64
    dtOld::Float64
    dt::Float64
    ab1::Float64; ab2::Float64
    cfl::Float64
    ts::String
    checkTimeStep::Int
    step::Int
    Nt::Int
    nuScaleFactor::Float64
    # --- BCs ---
    bcs::String
    bc::Matrix{Int}              # bc[side,axis]
    bcLabel::String
    multipleImplicitSolversNeeded::Bool
    # --- manufactured / known solution ---
    ms::String
    knownSolution::String
    degreex::Int; degreet::Int
    kx::Float64; ky::Float64; kt::Float64
    tzScale::Int
    computeErrors::Int
    # exact solution function handles (Any = unset until defined)
    ue::Any; ve::Any; pe::Any
    uex::Any; vex::Any; pex::Any
    uey::Any; vey::Any; pey::Any
    uet::Any; vet::Any; pet::Any
    uexx::Any; vexx::Any; pexx::Any
    uexy::Any; vexy::Any; pexy::Any
    ueyy::Any; veyy::Any; peyy::Any
    uett::Any; vett::Any; pett::Any
    ufe::Any; vfe::Any; pfe::Any
    # --- BC function handles ---
    guax::Any; gubx::Any; guay::Any; guby::Any
    gvax::Any; gvbx::Any; gvay::Any; gvby::Any
    pgax::Any; pgbx::Any; pgay::Any; pgby::Any
    gu::Vector{Any}   # {guax, gubx, guay, guby}
    gv::Vector{Any}
    gp::Vector{Any}
    u0::Any; v0::Any
    # --- initial conditions ---
    ic::String
    uic::Float64; vic::Float64
    shearBeta::Float64; shearDeltav::Float64
    # --- artificial dissipation ---
    ad::Int; ad21::Float64; ad22::Float64
    # --- outflow / inflow pressure ---
    outflowPressureCoeffp::Float64
    outflowPressureCoeffpn::Float64
    pOutflow::Float64
    pressureInflowValue::Float64
    uInflow::Float64
    # --- divergence damping ---
    cdv::Float64
    # --- sparse matrix factors ---
    dA::Any                  # LU for pressure
    dAimp::Vector{Any}       # LU(s) for implicit velocity
    # --- timing ---
    cpuSetup::Float64; cpuGetUt::Float64; cpuPressure::Float64
    cpuBC::Float64; cpuImplicit::Float64; cpuFactorImpMatrix::Float64
    cpuPlot::Float64; cpuTotal::Float64
    # --- output / diagnostics ---
    maxErr::Vector{Float64}
    uNorm::Vector{Float64}
    maxDivU::Float64; maxGradU::Float64
    # --- options ---
    idebug::Int
    plotOption::Int
    plotErrors::Int
    plotEveryStep::Int
    plotGrid::Int
    plotSolutionOnGhost::Int
    movieMode::Int
    savePlots::Int; savePlotThisStep::Bool
    figDir::String; plotName::String
    echo::Int
    checkFileName::String
    # --- curvilinear (annulus etc.) --- unused in Cartesian mode
    x0::Float64; y0::Float64
    startAngle::Float64; endAngle::Float64
    innerRadius::Float64; outerRadius::Float64
end

function Par()
    Par(
        # nd, nu, rho
        2, 0.1, 1.0,
        # xa,xb,ya,yb
        0.0,1.0, 0.0,1.0,
        # N0,Nx,Ny,Ngx,Ngy
        10,-1,-1, 0,0,
        # dx,dy
        0.0,0.0,
        # dr
        [0.0,0.0],
        # numGhost,orderInSpace
        1,2,
        # gid(2×2), dim(2×2)
        zeros(Int,2,2), zeros(Int,2,2),
        # isCartesian
        true,
        # x(0×0×2), rx(0×0×2×2)
        zeros(0,0,2), zeros(0,0,2,2),
        # map
        "Cartesian",
        # tf,tp,tFinal,dtMax,dtOld,dt,ab1,ab2,cfl,ts,checkTimeStep,step,Nt,nuScaleFactor
        0.5,0.1, 0.5, 1e10,-1.0,0.0,0.0,0.0, 0.9,"ab2",10000,1,0,1.0,
        # bcs,bc(2×2),bcLabel,multipleImplicitSolversNeeded
        "nnnn",zeros(Int,2,2),"",false,
        # ms,knownSolution,degreex,degreet,kx,ky,kt,tzScale,computeErrors
        "none","none",2,2, 1.0,1.0,0.5, 0, 0,
        # exact solution functions (unset)
        nothing,nothing,nothing,  # ue,ve,pe
        nothing,nothing,nothing,  # uex,vex,pex
        nothing,nothing,nothing,  # uey,vey,pey
        nothing,nothing,nothing,  # uet,vet,pet
        nothing,nothing,nothing,  # uexx,vexx,pexx
        nothing,nothing,nothing,  # uexy,vexy,pexy
        nothing,nothing,nothing,  # ueyy,veyy,peyy
        nothing,nothing,nothing,  # uett,vett,pett
        nothing,nothing,nothing,  # ufe,vfe,pfe
        # BC function handles (unset)
        nothing,nothing,nothing,nothing,  # guax,gubx,guay,guby
        nothing,nothing,nothing,nothing,  # gvax,gvbx,gvay,gvby
        nothing,nothing,nothing,nothing,  # pgax,pgbx,pgay,pgby
        Any[],Any[],Any[],  # gu,gv,gp
        nothing,nothing,  # u0,v0
        # ic,uic,vic,shearBeta,shearDeltav
        "default",1.0,0.0, 40.0,1e-2,
        # ad,ad21,ad22
        0,0.1,0.1,
        # outflowPressureCoeffp,outflowPressureCoeffpn,pOutflow,pressureInflowValue,uInflow
        1.0,1.0,0.0, 1.0,1.0,
        # cdv
        1.0,
        # dA,dAimp
        nothing, Any[],
        # timing
        0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,
        # maxErr,uNorm,maxDivU,maxGradU
        zeros(4), zeros(3), 0.0,0.0,
        # idebug,plotOption,plotErrors,plotEveryStep,plotGrid,plotSolutionOnGhost,movieMode
        0,1,1,0,0,0, 0,
        # savePlots,savePlotThisStep,figDir,plotName,echo
        0,false,"fig","ins",0,
        # checkFileName
        "ins.check",
        # curvilinear
        0.0,0.0, 0.0,1.0, 0.5,1.0
    )
end
