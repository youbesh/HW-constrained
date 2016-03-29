


module AssetTests

	using FactCheck, HW_constrained

	context("testing components") do

		facts("finite differences") do

			f(x,g) = 2.*sum(x)
			@fact HW_constrained.finite_diff(f,[4.0]) --> [2.0]
			@fact HW_constrained.finite_diff(f,[4.0;rand()]) --> [2.0;2.0]

			ff(x,g) = sum(x.^3)
			@fact HW_constrained.finite_diff(ff,[2.0])[1] --> roughly(12)
			@fact HW_constrained.finite_diff(ff,[2.0;1.0]) --> roughly([12.0;3.0])

		end

		facts("test_finite_diff") do
			function f(x,g)
				g[:] = 2
				return 2*sum(x)
			end
			@fact HW_constrained.test_finite_diff(f,[4.0],1e-4) --> true
			function f(x,g)
				g[:] = 5
				return 2*sum(x)
			end
			@fact HW_constrained.test_finite_diff(f,[4.0],1e-4) --> (false,[1])
			function f(x,g)
				g[:] = 5
				return 2*sum(x)
			end
			@fact HW_constrained.test_finite_diff(f,vcat(1.0,2,3,4,5),1e-4) --> (false,vcat(1,2,3,4,5))
		end


		facts("tests gradient of objective function") do
			d = HW_constrained.data(0.5)
			x = ones(Float64,d["na"]+1)
			@fact HW_constrained.test_finite_diff((ix,g)->HW_constrained.obj(ix,g,d),x) --> true

			x = rand(d["na"]+1)
			@fact HW_constrained.test_finite_diff((ix,g)->HW_constrained.obj(ix,g,d),x) --> true
		end


		facts("tests gradient of constraint function") do
			d = HW_constrained.data(0.5)
			x = ones(Float64,d["na"]+1)
			@fact HW_constrained.test_finite_diff((ix,g)->HW_constrained.constr(ix,g,d),x) --> true

			x = rand(d["na"]+1)
			@fact HW_constrained.test_finite_diff((ix,g)->HW_constrained.constr(ix,g,d),x) --> true
		end
	end

	context("testing result of both maximization methods") do

		truth = HW_constrained.DataFrame(a=[0.5;1.0;5.0],c = [1.008;1.004;1.0008],omega1=[-1.41237;-0.20618;0.758763],omega2=[0.801455;0.400728;0.0801455],omega3=[1.60291;0.801455;0.160291],fval=[-1.20821;-0.732819;-0.013422])

		facts("checking result of NLopt maximization") do

			t1 = table_NLopt()
			for c in names(truth)
				@fact t1[c] --> roughly(truth[c],atol=1.e-4)
			end
		end


		facts("checking result of NLopt maximization") do
			t1 = table_JuMP()
			for c in names(truth)
				@fact t1[c] --> roughly(truth[c],atol=1.e-4)
			end
		end
	end




end



