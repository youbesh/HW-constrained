

# constrained maximization exercises

## portfolio choice problem

module HW_constrained

	using JuMP, NLopt, DataFrames

	export data, table_NLopt, table_JuMP

	function data(a=0.5)
	
	end
	

	function max_JuMP(a=0.5)
	
	end

	function table_JuMP()
	
	end

	
	function obj(x::Vector,grad::Vector,data::Dict)
	
	end

	function constr(x::Vector,grad::Vector,data::Dict)
	  
	end

	function max_NLopt(a=0.5)
		
	end

	function table_NLopt()
	
	end

	# function `f` is for the NLopt interface, i.e.
	# it has 2 arguments `x` and `grad`, where `grad` is
	# modified in place
	# if you want to call `f` with more than those 2 args, you need to
	# specify an anonymous function as in
	# other_arg = 3.3
	# test_finite_diff((x,g)->f(x,g,other_arg), x )
	# this function cycles through all dimensions of `f` and applies
	# the finite differencing to each. it prints some nice output.
	function test_finite_diff(f::Function,x::Vector{Float64},tol=1e-6)
	
	end

	# do this for each dimension of x
	# low-level function doing the actual finite difference
	function finite_diff(f::Function,x::Vector)
		
	end

	function runAll()
		println("running tests:")
		include("test/runtests.jl")
		println("")
		println("JumP:")
		table_JuMP()
		println("")
		println("NLopt:")
		table_NLopt()
		ok = input("enter y to close this session.")
		if ok == "y"
			quit()
		end
	end


end


