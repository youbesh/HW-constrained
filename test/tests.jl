#LOAD_PATH
#push!(LOAD_PATH, "C:/Users/Youssef/Desktop/M2_semester_ 2/Numerical_Methods/HW-constrained/src")

module AssetTests

	using Base.Test, HW_constrained

	#@testset "testing components" begin

	#	@testset "finite differences" begin


	#	end

	#	@testset "test_finite_diff" begin
	#	end


	#	@testset "tests gradient of objective function" begin
	#	end


	#	@testset "tests gradient of constraint function" begin
	#	end
	#end

	@testset "testing result of both maximization methods" begin
	truth = DataFrame()
	truth[:a]=[0.5;1;5]
	truth[:c]=[1.00801;1.00401;1.0008]
	truth[:omega1]=[-1.41237;-0.206197;0.758762]
	truth[:omega2]=[0.801458;0.400729;0.0801456]
	truth[:omega3]=[1.60291; 0.801462; 0.160291]
	truth[:fvalue]=[-1.20821; -0.732819; -0.013422]

		@testset "checking result of NLopt maximization" begin
		@test isapprox(HW_constrained.table_NLopt()[:1],truth[:1], rtol=1e-3, atol=0)
		@test isapprox(HW_constrained.table_NLopt()[:2],truth[:2], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_NLopt()[:3],truth[:3], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_NLopt()[:4],truth[:4], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_NLopt()[:5],truth[:5], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_NLopt()[:6],truth[:6], rtol=1e-2, atol=0)
		end


		@testset "checking result of JUMP maximization" begin
		@test isapprox(HW_constrained.table_JuMP()[:1],truth[:1], rtol=1e-3, atol=0)
		@test isapprox(HW_constrained.table_JuMP()[:2],truth[:2], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_JuMP()[:3],truth[:3], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_JuMP()[:4],truth[:4], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_JuMP()[:5],truth[:5], rtol=1e-2, atol=0)
		@test isapprox(HW_constrained.table_JuMP()[:6],truth[:6], rtol=1e-2, atol=0)
		end
	end




end
