# constrained maximization exercises

## portfolio choice problem

module HW_constrained

	using JuMP, NLopt, DataFrames, Ipopt

	export data, table_NLopt, table_JuMP

	function dta(a)
	u(c) = -exp(-a * c)
	n = 3
	p = [1, 1, 1]
	e = [2, 0, 0]
	z2 = [0.72, 0.92, 1.12, 1.32]
	z3 = [0.86, 0.96, 1.06, 1.16]
	S=([0.72,0.86],[0.72,0.96],[0.72,1.06],[0.72,1.16],[0.92,0.86],[0.92,0.96],[0.92,1.06],[0.92,1.16],[1.12,0.86]
	,[1.12,0.96],[1.12,1.06],[1.12,1.16],[1.32,0.86],[1.32,0.96],[1.32,1.06],[1.32,1.16])
	dict = Dict("utility" => u, "a" => a, "numass" => n, "endw" => e,  "prices" => p,
	"return2" => z2, "return3" => z3, "states" => S)
	return dict
	end


	function obj(x::Vector,grad::Vector,dta::Dict)
		if length(grad) > 0
			grad[1] = dta["a"] * exp(-dta["a"] * x[1])

			grad[2] =  (1/16)*sum(
			dta["a"] * exp(
			-dta["a"] * (x[2]+x[3]*dta["states"][s][1] + x[4]*dta["states"][s][2])
			)
			for s = 1:16)

			grad[3] = (1/16)*sum(
			dta["a"] * exp(
			-dta["a"] * (x[2]+x[3]*dta["states"][s][1] + x[4]*dta["states"][s][2])
			) * dta["states"][s][1]
			for s = 1:16)

			grad[4] = (1/16) * sum(
			dta["a"] * exp(
			-dta["a"] * (x[2]+x[3]*dta["states"][s][1] + x[4]*dta["states"][s][2])
			) * dta["states"][s][2]
			for s = 1:16)
			end
			S = dta["states"]
			return -exp(-dta["a"]*x[1])+(1/16)*sum((-exp(-dta["a"]*(x[2]+x[3]*S[s][1]+x[4]*S[s][2]))) for s=1:16)
	end

	function constr(x::Vector,grad::Vector,dta::Dict)
		if length(grad) > 0
		grad[1] = 1
		grad[2] = 1
		grad[3] = 1
		grad[4] = 1
		end
		return x[1] + sum(x[i]- dta["endw"][i-1] for i in 2:4)
	end

	function max_NLopt(a)
		opt = Opt(:LD_MMA, 4)
		data = dta(a)
		#lower_bounds!(opt, [-Inf, 0.])
		xtol_rel!(opt,1e-4)
		max_objective!(opt, (x,grad) -> obj(x, grad, data))

		inequality_constraint!(opt, (x,grad) -> constr(x,grad,data))

		#ăcall optimize
		(optf,optx,ret) = optimize(opt, [1, -1.5, 1, 1.5])
	end

	function table_NLopt()
		df1=DataFrame()
		df1[:a]=[0.5;1;5]
		df1[:c]=[max_NLopt(0.5)[2][1];max_NLopt(1)[2][1];max_NLopt(5)[2][1]]
		df1[:omega1]=[max_NLopt(0.5)[2][2];max_NLopt(1)[2][2];max_NLopt(5)[2][2]]
		df1[:omega2]=[max_NLopt(0.5)[2][3];max_NLopt(1)[2][3];max_NLopt(5)[2][3]]
		df1[:omega3]=[max_NLopt(0.5)[2][4];max_NLopt(1)[2][4];max_NLopt(5)[2][4]]
		df1[:fvalue]=[max_NLopt(0.5)[1];max_NLopt(1)[1];max_NLopt(5)[1]]
	return df1
	end





	function max_JuMP(a)
		m=Model(solver = IpoptSolver())
		n=3
		@variable(m,c>=0)
		@variable(m,omega[1:n])
		pi=0.25*0.25
		S = dta(a)["states"]
		@NLobjective(m,Max,-exp(-a*c)+pi*sum(-exp(-a*(omega[1]+omega[2]*S[s][1]+omega[3]*S[s][2])) for s in 1:16))
		@NLconstraint(m, c+sum(omega[i]- dta(a)["endw"][i] for i in 1:n)==0)
		print(m)
		status=solve(m)
		v=[getobjectivevalue(m);getvalue(c);getvalue(omega)[1];getvalue(omega)[2];getvalue(omega)[3]]
		return v
	end

	function table_JuMP()
		df2=DataFrame()
		df2[:A]=[0.5;1;5]
		df2[:B]=[max_JuMP(0.5)[2];max_JuMP(1)[2];max_JuMP(5)[2]]
		df2[:C]=[max_JuMP(0.5)[3];max_JuMP(1)[3];max_JuMP(5)[3]]
		df2[:D]=[max_JuMP(0.5)[4];max_JuMP(1)[4];max_JuMP(5)[4]]
		df2[:E]=[max_JuMP(0.5)[5];max_JuMP(1)[5];max_JuMP(5)[5]]
		df2[:F]=[max_JuMP(0.5)[1];max_JuMP(1)[1];max_JuMP(5)[1]]
		newname2=["a","c","omega1","omega2","omega3","fvalue"]
		names!(df2.colindex, map(parse, newname2))
		return df2
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
	#function test_finite_diff(f::Function,x::Vector{Float64},tol=1e-6)

	#end

	# do this for each dimension of x
	# low-level function doing the actual finite difference
	#function finite_diff(f::Function,x::Vector)

	#end

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
