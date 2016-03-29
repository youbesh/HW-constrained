

# constrained maximization exercises

## portfolio choice problem

module HW_constrained

	using JuMP, NLopt, DataFrames

	export data, table_NLopt, table_JuMP

	function data(a=0.5)
		na = 3  # number of assets
		nc = na + 1  # number of choice variables: na assets + 1 consumption
		ns = 4 # number of states for each asset
		nss = ns^2 # number of total states of the world
		e = [2.0;0;0]
		p = [1.0;1;1]
		@assert length(e) == na
		z = hcat(ones(ns^2),repeat([0.72,0.92,1.12,1.32],inner=[ns],outer=[1]),repeat([0.86,0.96,1.06,1.16],inner=[1],outer=[ns]))
		pi = Float64[1/size(z,1) for i=1:size(z,1)]
		return Dict("a"=>a,"na"=>na,"nc"=>nc,"ns"=>ns,"nss"=>nss,"e"=>e,"p"=>p,"z"=>z,"pi"=>pi)
	end
	function data2(a=0.5)
		na = 3  # number of assets
		ns = 4 # number of states for each asset
		nss = ns^2 # number of total states of the world
		e = [2;0;0]
		p = [1;1;1]
		@assert length(e) == na
		z = hcat(ones(ns^2),repeat([0.72,0.92,1.12,1.32],inner=[ns],outer=[1]),repeat([0.86,0.96,1.06,1.16],inner=[1],outer=[ns]))
		pi = Float64[1/size(z,1) for i=1:size(z,1)]
		return (a,na,ns,nss,e,p,z,pi)
	end

	function max_JuMP(a=0.5)
		d = data(a)
		m = Model()
		@defVar(m,c >= 0.0)
		@defVar(m,omega[1:d["na"]])
		@setNLObjective(m,:Max,-exp(-a*c) + sum{d["pi"][s] * (-exp(-a * sum{omega[j]*d["z"][s,j],j=1:d["na"]})),s=1:d["nss"]})
		@addNLConstraint(m,c + sum{d["p"][i]*(omega[i]-d["e"][i]),i=1:d["na"]} == 0.0)
		print(m)
		status = solve(m)
		return Dict("obj"=>getObjectiveValue(m),"c"=>getValue(c),"omega"=>getValue(omega))
	end

	function table_JuMP()
		d = DataFrame(a=[0.5;1.0;5.0],c = zeros(3),omega1=zeros(3),omega2=zeros(3),omega3=zeros(3),fval=zeros(3))
		for i in 1:nrow(d)
			xx = max_JuMP(d[i,:a])
			d[i,:c] = xx["c"]
			d[i,:omega1] = xx["omega"][1]
			d[i,:omega2] = xx["omega"][2]
			d[i,:omega3] = xx["omega"][3]
			d[i,:fval] = xx["obj"]
		end
		return d
	end

	u(a_,c) = -exp(-a_ * c)
	up(a_,c) = a_*exp(-a_*c)
	# using Base.Test
	# @test -1.0 == u(0,12.1)
	# @test -1.0 == u(a,0)
	# @test_approx_eq -exp(1) u(a,-1/a)

	# function obj(omega::Vector,grad::Vector,a_,p_,e_,pi_::Vector,z_::Matrix)
	function obj(x::Vector,grad::Vector,data::Dict)
	    nz,na = size(data["z"])
		c = x[1]
		omega = x[2:end]
	    # println("obj: omega = $omega")
	    # println("obj: c= $c")
	    EV = 0.0
	    for is=1:nz	 # for each state
	        EV += data["pi"][is] * u(data["a"],dot(omega, vec(data["z"][is,:])))
	    end
	    if length(grad) > 0
	    	grad[1] = up(data["a"],c)
	    	for i=1:data["na"] 
		    	grad[i+1] = 0
		    	for is =1:nz
		    		# println("omega = $omega")
		    		# println("vec(")
		    		# println(vec(data["z"][is,:]))
		    		# println("dot = $(dot(omega,vec(data["z"][is,:])))")
		    		grad[i+1] += data["pi"][is] * up(data["a"],dot(omega,vec(data["z"][is,:]))) * data["z"][is,i]
		    	end
		    end
	    end    
	    out = u(data["a"],c) + EV
	    # println("objective = $out")
	    # println("obj: grad = $grad")
	    return out
	end

	function constr(x::Vector,grad::Vector,data::Dict)
	    nz,na = size(data["z"])
		c = x[1]
		omega = x[2:end]

		# needs to be of form mycontraint(x) \leq 0
	    constr_val = dot(data["p"], data["e"].-omega) - c

	    if length(grad) > 0
	        grad[1] = -1	# wrt c
	        grad[2:end] = -data["p"] 	# wrt omega
	    end
	    # println("constrain= $constr_val")
	    # println("constrain grad = $grad")
	    return constr_val 
	end

	function max_NLopt(a=0.5)
		d = data(a)
		opt = NLopt.Opt(:LD_SLSQP,d["nc"])
		max_objective!(opt,(x,g)->obj(x,g,d))
		lower_bounds!(opt,[0;[-Inf for i=1:d["na"]]])
		upper_bounds!(opt,[+Inf for i=1:d["nc"]])
		equality_constraint!(opt,(x,g)->constr(x,g,d),1e-5)
		ftol_rel!(opt,1e-9)
		(optf,optx,ret) = optimize(opt, rand(d["nc"]))
	end

	function table_NLopt()
		d = DataFrame(a=[0.5;1.0;5.0],c = zeros(3),omega1=zeros(3),omega2=zeros(3),omega3=zeros(3),fval=zeros(3))
		for i in 1:nrow(d)
			xx = max_NLopt(d[i,:a])
			for j in 2:ncol(d)-1
				d[i,j] = xx[2][j-1]
			end
			d[i,end] = xx[1]
		end
		return d
	end

	function test_finite_diff(f::Function,x::Vector{Float64},tol=1e-6)
		# get gradient from f
		grad = similar(x)
		y = f(x,grad)
		# get finite difference approx
		fdiff = finite_diff(f,x)
		r = hcat(1:length(x),grad,fdiff,abs(grad-fdiff))
		errors = find(abs(grad-fdiff).>tol)
		if length(errors) >0
			println("elements with errors:")
			println("id  supplied gradient     finite difference     abs diff")
			for i in 1:length(errors)
				@printf("%d   %f3.8            %f3.8          %f1.8\n",r[errors[i],1],r[i,2],r[i,3],r[i,4])
			end
			return (false,errors)
		else 
			println("no errors.")
			return true
		end
	end

	# do this for each dimension of x
	function finite_diff(f::Function,x::Vector)
		h = sqrt(eps())
		fgrad = similar(x)
		tgrad = similar(x)
		for i in 1:length(x)
			step = abs(x[i]) > 1 ? abs(x[i]) * h : 1.0 * h
			newx = copy(x)
			newx[i] = x[i]+step
			fgrad[i] = (f(newx,tgrad) - f(x,tgrad))/step
		end
		return fgrad
	end

	"""
	Version for two separate functions supplying objective and gradient
	"""
	function test_finite_diff(f::Function,g::Function,x::Vector{Float64},tol=1e-6)
		# get gradient from f
		grad = g(x)
		# get finite difference approx
		fdiff = finite_diff2(f,x)
		r = hcat(1:length(x),grad,fdiff,abs(grad-fdiff))
		errors = find(abs(grad-fdiff).>tol)
		if length(errors) >0
			println("elements with errors:")
			println("id  supplied gradient     finite difference     abs diff")
			for i in 1:length(errors)
				@printf("%d   %f3.8            %f3.8          %f1.8\n",r[errors[i],1],r[i,2],r[i,3],r[i,4])
			end
			return (false,errors)
		else 
			println("no errors.")
			return true
		end
	end
	function finite_diff2(f::Function,x::Vector)
		h = sqrt(eps())
		fgrad = similar(x)
		tgrad = similar(x)
		for i in 1:length(x)
			step = abs(x[i]) > 1 ? abs(x[i]) * h : 1.0 * h
			newx = copy(x)
			newx[i] = x[i]+step
			fgrad[i] = (f(newx) - f(x))/step
		end
		return fgrad
	end


	function runAll()
		println("running tests:")
		include("test/runtests.jl")
		println("")
		ok = input("enter y to close this session.")
		if ok == "y"
			quit()
		end
	end


end


