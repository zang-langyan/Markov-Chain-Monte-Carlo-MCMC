using Distributions
using Random

""" # Configuration for MCMC
    MCMC_Config([chain::Int, init::Real, jump, posterior, space::AbstractVector, burnin::Int, rng]) -> MCMC_Config
* `chain::Int` - the length of Markov Chain to be generated (default 10000)
* `init::Real` - initial value for the chain of θ (default 0.5)
* `jump` - the distribution that proposed jump (Δθ) follows (default Distributions.Normal(0,0.2))
* `posterior` 
    - (`::Symbol`) the Symbol of self-defined density function
    - (`::var(#)`) a anonymous function of the density function (default θ -> Distributions.pdf(Beta(15,7),θ))
* `space::AbstractVector` - the accepted span of the parameter θ [min,max] (default [0,1])
* `burnin::Int` - the length of the chain to be burned from the beginning (default 0)
* `rng` - random seed (default 42), set `rng = nothing` to use no random seed
"""
Base.@kwdef struct MCMC_Config
    chain::Int = 10000
    init::Real = 0.5
    jump = Distributions.Normal(0,0.2)
    posterior = θ -> Distributions.pdf(Beta(15,7),θ) # posterior = :p
    space::AbstractVector = [0,1]
    burnin::Int = 0
    rng = 42
end
""" # Metropolis Algorithm
    metropolis(dfunc, chain::Int, theta_init::Real; [jump = Normal(0,0.2), space_min::Real = -Inf, space_max::Real = Inf, burnin::Int = 0, rng = nothing]) -> AbstractVector

Compute a Markov Chain using Metropolis Algorithm

# Parameters
## arguments
* `dfunc`
    - (`::Symbol`) the Symbol of self-defined density function
    - (`::var(#)`) a anonymous function of the density function
* `chain::Int` - the length of Markov Chain to be generated
* `theta_init::Real` - initial value for the chain of θ (default 0.5)

## keyword arguments
* `jump` - the distribution that proposed jump (Δθ) follows (default Distributions.Normal(0,0.2))
* `space_min::Real` - the accepted minimum of the parameter θ (default -Inf)
* `space_max::Real` - the accepted maximum of the parameter θ (default  Inf)
* `burnin::Int` - the length of the chain to be burned from the beginning (default 0)
* `rng` - random seed (default nothing), set `rng = 0` to use seed `MersenneTwister(0)`

# Examples
```julia-repo
julia> p(θ) = Distributions.pdf(Beta(15,7),θ) # define a density function

julia> metropolis(:p, 5000, 0.5; jump = Normal(0,0.2), space_min = 0, space_max = 1, burnin = 1000, rng = 42)
4000-element Vector{Float64}:
 0.7254870615709366
 0.5374079418869023
 0.5374079418869023
 0.5161750364065404
 0.7377560837277216
 0.7377560837277216
 0.8127557482035582
 ⋮
 0.7166820320569407
 0.7166820320569407
 0.7166820320569407
 0.7166820320569407
 0.7166820320569407
 0.7586353644413187

julia> @time metropolis(θ -> Distributions.pdf(Beta(15,7),θ), 10000, 0.5; jump = Normal(0,0.5), space_min = 0, space_max = 1, burnin = 1000, rng = nothing)
0.061715 seconds (103.38 k allocations: 3.814 MiB, 88.90% compilation time)
9000-element Vector{Float64}:
 0.6057362032299476
 0.6057362032299476
 0.6057362032299476
 0.6057362032299476
 0.6057362032299476
 0.6057362032299476
 0.6057362032299476
 ⋮
 0.4516105573548591
 0.4516105573548591
 0.4516105573548591
 0.5014099540247418
 0.6720304713998495
 0.6720304713998495

```
"""
function metropolis(dfunc, chain::Int, theta_init::Real; jump = Normal(0,0.2), space_min::Real = -Inf, space_max::Real = Inf, burnin::Int = 0, rng = nothing)
    posterior = 0 # initialize posterior to the outer scope
    try
        if typeof(dfunc) == Symbol
            posterior = getfield(Main,dfunc)
        else
            posterior = dfunc
        end
    catch
        print("dfunc must be a self-defined function name or an anonymous")
    end

    θ_cur = theta_init
    θ_freq = [theta_init]

    Random.seed!(MersenneTwister(rng))

    while true
        Δθ = rand(jump)
        θ_pro = θ_cur + Δθ

        if θ_pro < space_min || θ_pro > space_max
            pmoving = 0
        elseif posterior(θ_cur) == 0
            pmoving = 1
        else
            pmoving = min(1, posterior(θ_pro)/posterior(θ_cur))
        end

        if rand() <= pmoving
            push!(θ_freq, θ_pro)
            θ_cur = θ_pro
        else
            push!(θ_freq, θ_cur)
        end

        if length(θ_freq) >= chain
            break
        end
    end

    return θ_freq[burnin + 1 : end]
end

""" 
    metropolis(Config::MCMC_Config) -> AbstractVector

Compute a Markov Chain using Metropolis Algorithm

use the struct constructor `::MCMC_Config` to configure the parameter for Metropolis Algorithm

check struct `MCMC_Config` for more details

# Examples
```julia-repo
julia> Config = MCMC_Config()

julia> metropolis(Config)
10000-element Vector{Float64}:
 0.5
 0.6828984833360351
 0.6828984833360351
 0.6655360165668127
 0.6655360165668127
 0.6655360165668127
 0.6591963728158979
 ⋮
 0.7246337711515528
 0.5807104209871898
 0.5275540714270037
 0.5569823211525173
 0.5272108581003658
 0.5272108581003658

julia> p(θ) = Distributions.pdf(Gamma(1,2),θ) # define a density function

julia> Config = MCMC_Config(chain = 5000, jump = Distributions.TDist(5), space = [0,Inf], posterior = :p, burnin = 500, rng = nothing)

julia> metropolis(Config)
4500-element Vector{Float64}:
 0.34362723284778773
 0.34362723284778773
 1.1023030285120496
 1.1023030285120496
 1.1023030285120496
 1.1023030285120496
 1.565816384278301
 ⋮
 2.05294789570401
 0.80709730394856
 0.80709730394856
 0.9931268934687827
 0.9931268934687827
 0.9931268934687827

```
"""
function metropolis(Config::MCMC_Config)
    dfunc = Config.posterior
    chain = Config.chain
    theta_init = Config.init
    jump = Config.jump
    space_min = Config.space[1]
    space_max = Config.space[2]
    burnin = Config.burnin
    rng = Config.rng
    metropolis(dfunc, chain, theta_init; jump = jump, space_min = space_min, space_max = space_max, burnin = burnin, rng = rng)
end