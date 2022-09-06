"""
Module contains the step (aka go function) for the model
"""

module go
    using Agents
    using Random
    using StatsBase
    using DataFrames
    using Distributions

    include("Disturbance_functions.jl")
    
    #%This is the step function for the individual trees (no globals changed) 
    #* RUN ONCE PER AGENT PER CALL (I.E. Multiple times per tick if multiple agents)
    function agent_step!(
        agent,
        model
    )
        #agent.age += 1

        #disturbance_functions.lsp_disturbance(model)
    end

    #%This is the step function for global level changes e.g. ticks
    #* RUN ONCE PER MODEL PER CALL (I.E. ONCE PER TICK)
    function model_step!(model)
        model.tick += 1

        if model.disturbance_freq > 0
            disturbance_functions.lsp_disturbance(model)
        end
        
    end
end