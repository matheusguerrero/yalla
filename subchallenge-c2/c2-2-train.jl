cd("subchallenge-c2")

using NeuralEstimators
using Flux
using CSV
using DataFrames
using Tables
using Random: seed!
using RData


print(ARGS)
K=350000 # Traing data sixe
input=load(string("NBE/replicates/train_N",K,".Rdata")) 

# Normalise the data. Makes training easier

m=29.176
s=22.141

for i ∈ eachindex(input["Z_train"])
		input["Z_train"][i]=reshape(input["Z_train"][i],1,21000)
		input["Z_train"][i]=(input["Z_train"][i] .- m) ./ s

end
for i ∈ eachindex(input["Z_test"])
		input["Z_test"][i]=reshape(input["Z_test"][i],1,21000)
		input["Z_test"][i]=(input["Z_test"][i] .- m) ./ s
end
for i ∈ eachindex(input["Z_val"])
		input["Z_val"][i]=reshape(input["Z_val"][i],1,21000)
		input["Z_val"][i]=(input["Z_val"][i] .- m) ./ s

end

# ---- Architecture ---

seed!(1)
using Distributions


ψ = Chain(
	Dense(1, 48 ),
	Dense(48, 48 , relu)

)


ϕ = Chain(
	Dense(w1,w2,relu),
	Dense(w2, p)
)


print(ϕ)
θ̂ = DeepSet(ψ, ϕ)

# ---- Training ----
 

int_path = "NBE/intermediates/"
if !isdir(int_path) mkpath(int_path) end
savepath = "$int_path/runs"


# Define custome loss function

using Flux: mean
using Flux: Zygote
using ..Flux: ofeltype, epseltype
function _check_sizes(ŷ::AbstractArray, y::AbstractArray)
  for d in 1:max(ndims(ŷ), ndims(y)) 
   size(ŷ,d) == size(y,d) || throw(DimensionMismatch(
      "loss function expects size(ŷ) = $(size(ŷ)) to match size(y) = $(size(y))"
    ))
  end
end

function custom_loss(ŷ, y; agg = mean)
   #_check_sizes(ŷ, y)
   error1 = y .* ofeltype(ŷ, 0.99)  .- ŷ
   temp1 = Zygote.ignore_derivatives(error1 .>  0)
   error2 = ŷ.- y .* ofeltype(ŷ, 1.01) 
   temp2 = Zygote.ignore_derivatives(error2 .>  0)
   
   agg((ofeltype(ŷ, 0.9).*error1 ) .* temp1 .+  (ofeltype(ŷ, 0.1).*error2 ) .* temp2 )
end


using Flux: loadparams!


# Train the estimator
estimator = train(θ̂, input["theta_train"], input["theta_val"], input["Z_train"], input["Z_val"],savepath = savepath,loss=custom_loss,stopping_epochs=10, batchsize=64)






using Flux: loadparams!
loadparams!(θ̂, loadbestweights(string(savepath)))
loadparams!(estimator, loadbestweights(string(savepath)))

# ---- Assess the estimator ----


assessment = assess([estimator], input["theta_test"],input["Z_test"])
test_risk=risk(assessment,loss=custom_loss)
print(test_risk)

CSV.write(savepath * "/estimates.csv", assessment.df)
CSV.write(savepath * "/runtime.csv", assessment.runtime)#

input=CSV.read("../Data/Amaurot.csv", DataFrame)

tmp=input[1:end,"Y"]
tmp=(reshape(tmp,1,21000) .- m) ./ s

print("Estimate from data")
print(θ̂(tmp))


boot = bootstrap(θ̂, tmp; B = 1000)
savepath = "$int_path/estimates"
if !isdir(savepath) mkpath(savepath) end

import Tables: table

CSV.write(savepath * "/boot_estimates.csv", Tables.table(boot))



