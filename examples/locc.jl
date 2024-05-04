using Yao

# Follow Nielson & Chuang 12.5.1
# Alice do measurement: M1 = cos(θ) |0⟩⟨0| + sin(θ) |1⟩⟨1|, M2 = sin(θ) |0⟩⟨0| + cos(θ) |1⟩⟨1|
# where θ = atan(0.6/0.8) =  0.6435011087932843...
# Alice sends measurement result to Bob.
# They will both apply X gate if Alice's result is 2.
# Following is code to verify the result.

bell_state = ghz_state(2)
statevec(bell_state)
target_state = 0.8 * product_state(bit"00") + 0.6 * product_state(bit"11")

θ = atan(0.6 / 0.8)
M1 = ComplexF64[cos(θ) 0.0; 0.0 sin(θ)]
M2 = ComplexF64[sin(θ) 0.0; 0.0 cos(θ)]

meas1 = GeneralMatrixBlock(M1)
meas2 = GeneralMatrixBlock(M2)

res_state1 = copy(bell_state) |> kron(meas1, I2) |> normalize!
@assert res_state1 ≈ target_state
print_table(res_state1)
# Printed result
# 00 ₍₂₎   0.8 + 0.0im
# 01 ₍₂₎   0.0 + 0.0im
# 10 ₍₂₎   0.0 + 0.0im
# 11 ₍₂₎   0.6 + 0.0im

res_state2 = copy(bell_state) |> kron(meas2, I2) |> kron(X, X) |> normalize!
@assert res_state2 ≈ target_state
print_table(res_state2)
# Printed result
# 00 ₍₂₎   0.8 + 0.0im
# 01 ₍₂₎   0.0 + 0.0im
# 10 ₍₂₎   0.0 + 0.0im
# 11 ₍₂₎   0.6 + 0.0im