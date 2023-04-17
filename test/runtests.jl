using Test
using OrdinaryDiffEq: ODEProblem
using CompartmentalModelServer: greet, problem, solution, simulate, ModelSolution

@testset "CompartmentalModelServer" begin
    @testset "greet" begin @test greet() == "Hello World!" end
    @testset "problem" begin @test problem() isa ODEProblem end
    @testset "solution" begin @test solution() isa ModelSolution end
    @testset "simulate" begin @test simulate(; Snew = 888.0::Float64,
                                             βnew = 0.7::Float64) isa
                                    ModelSolution end
end
