using Test
using OrdinaryDiffEq
using CompartmentalModelServer: greet, problem, solution, simulate

@testset "CompartmentalModelServer" begin
    @testset "greet" begin
        @test greet() == "Hello World!"
    end
    @testset "problem" begin
        prob = problem()
        @test prob isa ODEProblem
    end
    @testset "solution" begin
        sol = solution()
        @test sol isa CompartmentalModelServer.ModelSolution
    end
    @testset "simulate" begin
        newsim = simulate(; Snew = 888.0::Float64, βnew = 0.7::Float64)
        @test newsim isa CompartmentalModelServer.ModelSolution
    end
end
