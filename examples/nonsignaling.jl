# %%
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra

# %% [markdown]
# We can represent the 4 possible functions with 1 bit input and 1 bit output as the following classical channels:
# assume f: x → y where x, y ∈ {0, 1}
# let us represent the input and output by the density matrix of the quantum state
# the classical channel is then a mapping from density matrix to density matrix
# i.e f1: 0/1 → 0 maps $\ket{0}\bra{0}$ to $\ket{0}\bra{0}$ and $\ket{1}\bra{1}$ to $\ket{0}\bra{0}$
# i.e f2: 0/1 → 1 maps $\ket{0}\bra{0}$ to $\ket{1}\bra{1}$ and $\ket{1}\bra{1}$ to $\ket{1}\bra{1}$
# i.e f3: 0/1 → 0/1 maps $\ket{0}\bra{0}$ to $\ket{0}\bra{0}$ and $\ket{1}\bra{1}$ to $\ket{1}\bra{1}$
# i.e f4: 0/1 → 1/0 maps $\ket{0}\bra{0}$ to $\ket{1}\bra{1}$ and $\ket{1}\bra{1}$ to $\ket{0}\bra{0}$
# we could represent the classical channel with a 4x4 matrix where the interpretation is as follows:
# entry of f_i to represent the mapping from $|x_i\rangle\langle x_i|$ to $|y_o\rangle\langle y_o|$ is given  by $|x_iy_o\rangle\langle x_iy_o|$ in the matrix

# %% 
f1 = [1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0]
f2 = [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1]
f3 = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]
f4 = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]

# %% [markdown]
# We denote whether a function is balanced or not by a single qubit. In the quantum circuit that implements the algorithm to test if a function is balanced or not, the output qubit will be in state $|0\rangle\langle 0|$ if the function is constant and $|1\rangle\langle 1|$ if the function is balanced. Hence, the target function we try to maximize is

# %%
function target(Mis)
	s0 = [1 0; 0 0]
	s1 = [0 0; 0 1]
	return (tr(s0 * Mis[1]*s0) + tr(s1 * Mis[2]*s0) + tr(s1 * Mis[3]*s0) + tr(s0 * Mis[4]*s0)) /4.0
end

J_pi = ComplexVariable(16,16)
J_2 = ComplexVariable(16,16)

kron(J_pi,[1 0; 0 1])
partialtrace(partialtrace(J_pi, 1, [2,2,2,2]),1,[2,2,2])

function deutsch(n_a::Int, n_b::Int,optimizer = Mosek.Optimizer,silent=true, signaling = false)

	J_pi = ComplexVariable(2^n_a *2^n_b,2^n_a*2^n_b)
	
	reshape(J_pi,2,2,2,2)

end


function nonsignaling_game(optimizer = Mosek.Optimizer, silent = true, signaling = false)
	# optimizer = SCS.Optimizer
	model = Model(optimizer)
	# silent = true
	silent && set_silent(model)

	pabxy = @variable(model, [1:2, 1:2, 1:2, 1:2], lower_bound = 0.0, upper_bound = 1.0)
	pax = @variable(model, [1:2, 1:2], lower_bound = 0.0, upper_bound = 1.0)
	pby = @variable(model, [1:2, 1:2], lower_bound = 0.0, upper_bound = 1.0)

	for given_val in 1:2
		@constraint(model, sum(pax[:, given_val]) == 1.0)
		@constraint(model, sum(pby[:, given_val]) == 1.0)
	end

	for y in 1:2
		for x in 1:2
			@constraint(model, sum(pabxy[:, :, x, y]) == 1.0)
			@constraint(model, sum(pabxy[:, :, x, y]) == 1.0)
			if !signaling
				for a in 1:2
					@constraint(model, sum(pabxy[a, :, x, y]) == pax[a, x])
				end
				for b in 1:2
					@constraint(model, sum(pabxy[:, b, x, y]) == pby[b, y])
				end
			end
		end
	end

	@objective(model, Max, sum([cabxy(a, b, x, y) * pabxy[a + 1, b + 1, x + 1, y + 1] for a in 0:1, b in 0:1, x in 0:1, y in 0:1]))

	optimize!(model)

	@assert is_solved_and_feasible(model)

	println(solution_summary(model))

    return value.(pabxy)
end