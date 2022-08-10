## TODO Write a description

module distribution_functions
    using Distributions
    ## TODO Document
    function generate_LogNormal(m,std)
        γ = 1+std^2/m^2
        μ = log(m/sqrt(γ))
        σ = sqrt(log(γ))

        return LogNormal(μ,σ)
    end
end