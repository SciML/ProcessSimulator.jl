using Clapeyron

#= n_feed = [5.0, 20.0, 5.0, 5.0, 5.0]
ν = [4  0; -5 -1; 1 -3; -3  1; 1  1]
ξ1 = -1.0:0.5:1.0
ξ2 = -4.0:2.0:1.0 
xgrid = ξ1'.*ones(size(ξ2, 1)) |> vec
ygrid =  ones(size(ξ1, 1))'.*ξ2 |> vec
ξ_samples = [xgrid ygrid]
n_samples = n_feed .+ ν*ξ_samples[1:end, :]'
is_negative = any(n_samples .< 0, dims = 1)
n_samples = n_samples[:, findall(x -> x == 0, is_negative[1:end])] =#


function Gibbs(T, P, n)
    comps =  ["carbon dioxide", "carbon monoxide", "hydrogen", "water", "methane"]
    model = ReidIdeal(comps;
    userlocations = (a = [19.8, 29443.7207031e-3, 23052.640625e-3, 33444.6210938e-3, 37.981], 
                    b = [73.4400024414e-3, -5.6729388237e-3, 33.7491416931e-3, -5.79920578003e-3, -74.6220016479e-3], 
                    c = [-0.0560199990869e-3, 0.0158794000745e-3, -0.0639906972647e-3, 0.025168100372e-3, 0.30189999938e-3],
                    d = [1.71499996213e-08,-6.43621388008e-09, 5.10229983774e-08, -1.43102997754e-08, -0.000283269997453e-3],
                    e = [0.0, 0.0, -1.3769899887e-11, 2.76249001452e-12, 9.07109978243e-11])) #Same used in COCO
   
    ∑n = sum(n)

    return gibbs_free_energy(model, P, T, n)/∑n, entropy(model, P, T, n)/∑n

end

Gibbs(500, 101325, [5.0, 20.0, 5.0, 5.0, 5.0])[2] - Gibbs(500, 101325, [1.0, 29.0, 16.0, 4.0, 0.0])[2]

-7595.1393 - (-5691.8746)

28.455021 - 23.447983

