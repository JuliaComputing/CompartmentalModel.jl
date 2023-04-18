module CompartmentalModelServer

using ModelingToolkit
using OrdinaryDiffEq
using StructTypes

export greet
export problem
export solution
export simulate

greet() = "Hello World!"

"Model default values of type `NamedTuple{Float64}`."
const DEFAULTS = (; S = 999.0::Float64, I = 1.0::Float64, R = 0.0::Float64,
                  β = (1 / 2)::Float64, γ = (1 / 3)::Float64,
                  tstart = 0.0::Float64, tstop = 30.0::Float64)

"Construct the compartmental model."
function system()::ODESystem

    # define independent variable and differentiable
    @variables t
    D = Differential(t)

    # define states and parameters
    states = @variables begin
        S(t) = DEFAULTS.S
        I(t) = DEFAULTS.I
        R(t) = DEFAULTS.R
    end
    params = @parameters begin
        β = DEFAULTS.β
        γ = DEFAULTS.γ
    end

    # define equations
    eqns = Equation[D(S) ~ -β * I * S / (S + I + R),
                    D(I) ~ (β * I * S / (S + I + R)) - (γ * I),
                    D(R) ~ γ * I]

    # construct and return ode system
    @named _model = ODESystem(eqns, t, states, params)
    return structural_simplify(_model)
end

"Constructs the original compartmental model ode problem."
function problem(; states = Float64[], params = Float64[],
                 tspan = (DEFAULTS.tstart, DEFAULTS.tstop))
    ODEProblem(system(), states, tspan, params)
end

"Holds values for beta and gamma as well as serialization/deserialization information."
struct ModelParameters
    beta::Float64
    gamma::Float64
end

"Holds values complete model solution as well as serialization/deserialization information."
struct ModelSolution
    time::Vector{Float64}
    S::Vector{Float64}
    I::Vector{Float64}
    R::Vector{Float64}
    params::ModelParameters
end

function ModelSolution(sol::ODESolution, prob::ODEProblem)
    ModelSolution(sol.t,
                  sol[:S], sol[:I], sol[:R],
                  ModelParameters(prob.p...))
end

StructTypes.StructType(::Type{ModelSolution}) = StructTypes.Struct()

"Constructs the original compartmental model ode solution and returns an object that can be serialized to JSON."
solution = (prob = problem()) -> ModelSolution(solve(prob, Tsit5()), prob)

"""
    simulate(; Snew, Inew, Rnew, βnew, γnew, tstart, tstop)::ModelSolution

Returns a new solution for the compartmental model given the initial conditions, parameters
and timespan passed as keyword arguments. The result is a `ModelSolution` object that can be
readily serialized to JSON.
"""
function simulate(; Snew = DEFAULTS.S, Inew = DEFAULTS.I, Rnew = DEFAULTS.R,
                  βnew = DEFAULTS.β, γnew = DEFAULTS.γ,
                  tstart = DEFAULTS.tstart, tstop = DEFAULTS.tstop)

    # use originally defined model
    sys = system()

    # unpack states and parameters
    t = ModelingToolkit.get_iv(sys)
    @unpack S(t), I(t), R(t) = ModelingToolkit.get_states(sys)
    β, γ = ModelingToolkit.get_ps(sys)

    # construct new problem
    newprob = remake(problem(); u0 = [Snew, Inew, Rnew], tspan = (tstart, tstop),
                     p = [β => βnew, γ => γnew])

    # solve and return
    newsol = solve(newprob, Tsit5())
    return ModelSolution(newsol, newprob)
end

function __init__()
    system()
    problem()
    solution()
    simulate()
end

end # module CompartmentalModelServer
