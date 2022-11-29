"""
Module contains all functions related to the spread and impact of the different disease types
implemented in the model.
"""

module disease_functions
    using Agents

    #//-------------------------------------------------------------------------------------------#
    #% PHYTOTHERA SPREAD
    """
    Function implements the spread of a soil-borne Phytothera like disease, e.g. Kauri dieback. The
    spread function follows a similar design to a number of implemented Phytothera (and other 
    disease) spread models, e.g. Harwood et al. (2009). This disease spread model is broadly a
    SEI (Susceptible, Exposed, Infected) model, with no recovery or immunity. The function contains
    two main exposure methods, global and local both of which may only act on a susceptible target
    species, defined in the demography input file. The global exposure method is a simple random
    exposure, where the probability of exposure is defined by the `phyto_global_infection_prob`
    user input. The local exposure probability is density dependent, with each infected neighbour
    tree within the `phyto_infectious_radius` having an equal chance of infecting a susceptible
    target tree. The probability of infection is defined by the `phyto_local_infection_prob` user
    input.

    ## Arguments
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the 
    current model being exposed to the disease.
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `global_infection_prob::Float64`: Probability of a susceptible target tree being infected by 
    the disease, if exposed to the disease by the global exposure method.
    - `phyto_local_prob::Float64`: Probability of a susceptible target tree being infected by the
    disease, if exposed by one of its infected neighbours. Note this probability is retested
    for each infected neighbour.
    - `infectious_radius::Int64`: Radius of the local exposure method, i.e. the number of
    cells around the target tree that may infect it.
    - `transmission_age::Int64`: Age (i.e. ticks since infection) at which a tree is capable of 
    transimitting the disease to the target tree. This is tested for each infected neighbour.
    - `pos::Tuple{Int64,Int64}`: Position of the target tree in the model grid.
    """
    function phytothera_spread(agent,
                              model,
                              global_infection_prob::Float64,
                              phyto_local_prob::Float64,
                              infectious_radius::Int64,
                              transmission_age::Int64,
                              pos::Tuple{Int64, Int64},)
        #* Check whether agent gets infected by global chance
        #TODO add an increased chance if the agent is an edge patch
        if rand() < global_infection_prob
            agent.phytothera_infected::Bool = true
            
        #* If not globally infected check whether agent gets infected by local neighbours 
        else
            for i in collect(nearby_ids(pos, model::ABM, infectious_radius))
                if agent.phytothera_infected == true
                    break
                else model[i].phytothera_infected_age ≥ transmission_age && rand() < phyto_local_prob
                        agent.phytothera_infected::Bool = true
                end
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% PHYTOTHERA IMPACT
    """
    This function checks to see whether a tree infected with a Phytothera like disease has become
    symptomatic and if so tests whether the tree dies or not. The probility of a tree becoming
    symptomatic is defined by the `phyto_symptoms_dev_prob` user input. The probability of a tree
    dying once symptomatic is defined by the `phyto_mortality_prob` user input.

    ## Arguments
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the
    current model which has been infected by the disease.
    - `model`: The AgentBasedModel object defining the current model. This object is usually
    created using `Agents.ABM()`.
    - `phytothera_infected::Bool`: Boolean value indicating whether the tree has been infected
    by the disease.
    - `min_symptomic_age::Int64`: Minimum age (i.e. ticks since infection) at which a tree may
    become symptomatic.
    - `symptom_prob::Float64`: Probability of a tree becoming symptomatic once infected and passing
    the minimum symptomatic age.
    - `mortality_prob::Float64`: Probability of a tree dying each tick once symptomatic.
    - `gap_maker::Int64`: Species specific property indicating whether a species is capable 
    of creating a forest gap (value = 1) or not (value = 0).
    - `species_ID::Int64`: Selected species ID
    - `cell::Int64`: Patch ID of cell where target tree is located
    - `expand::BitVector`: Flag indicating that a patch should be checked for gap expansion 
    (value = 1) or not (value = 0). See `expand_gap()`
    - `previous_species::Vector{Float64}`: Species ID of the last tree to have occupied every cell 
    in the current model grid.
    - `previous_height::Vector{Float64}`: Final height of the last tree to have occupied every cell 
    in the current model grid.
    - `a_height::Float64`: Height in meters of the tree.
    - `id::Int64`: Agent ID of the selected tree/agent
    """
    function phytothera_impact(agent,
                               model,
                               phytothera_infected::Bool,
                               min_symptomic_age::Int64,
                               symptom_prob::Float64,
                               mortality_prob::Float64,
                               gap_maker::Int64,
                               species_ID::Int64,
                               cell::Int64,
                               expand::BitVector,
                               previous_species::Vector{Float64},
                               previous_height::Vector{Float64},
                               a_height::Float64,
                               id::Int64,
                               )

        if phytothera_infected == true
            agent.phytothera_infected_age += 1
        end

        #* Check whether tree has become symptomatic
        if agent.phytothera_infected_age ≥ min_symptomic_age && 
            agent.phytothera_symptomatic == false
            if rand() < symptom_prob
                agent.phytothera_symptomatic::Bool = true
            end
        end

        #* If symptomatic chance for tree to die
        if agent.phytothera_symptomatic == true 
            if rand() < mortality_prob
                if gap_maker[species_ID] == 1
                    expand[cell] = true
                end
    
                ## Record dying tree species and height as a cell list
                previous_species[cell] = species_ID
                previous_height[cell] = a_height
    
                kill_agent!(id, model)
            end    
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% RUST FUNGUS SPREAD
    """
    Function implements the spread of an air-borne pathogen such as a rust e.g. Myrtle rust. The
    spread function follows a similar, simplified, design as that implemented for Phytothera. The 
    function contains only a exposure method where the probability of exposure is defined by the 
    rust_global_infection_prob`user input. 

    ## Arguments
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the 
    current model being exposed to the disease.
    - `rust_global_infection_prob`: Probability of a susceptible target tree being infected by 
    the disease, if exposed to the disease by the global exposure method.
    """
    function rust_spread(agent,
                         rust_infection_prob::Float64)
        #* Check whether agent gets infected by global chance
        #TODO add an increased chance if the agent is an edge patch
        if rand() < rust_infection_prob
            agent.rust_infected::Bool = true
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% RUST FUNGUS IMPACT
    """
    This function checks to see whether a tree infected with a rust like disease has become
    symptomatic and if so tests whether the tree dies or not. The probility of a tree becoming
    symptomatic is defined by the `rust_symptoms_dev_prob` user input. The probability of a tree
    dying once symptomatic is defined by the `rust_mortality_prob` user input.

    ## Arguments
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the
    current model which has been infected by the disease.
    - `model`: The AgentBasedModel object defining the current model. This object is usually
    created using `Agents.ABM()`.
    - `rust_infected::Bool`: Boolean value indicating whether the tree has been infected
    by the disease.
    - `min_symptomic_age::Int64`: Minimum age (i.e. ticks since infection) at which a tree may
    become symptomatic.
    - `symptom_prob::Float64`: Probability of a tree becoming symptomatic once infected and passing
    the minimum symptomatic age.
    - `mortality_prob::Float64`: Probability of a tree dying each tick once symptomatic.
    - `gap_maker::Int64`: Species specific property indicating whether a species is capable 
    of creating a forest gap (value = 1) or not (value = 0).
    - `species_ID::Int64`: Selected species ID
    - `cell::Int64`: Patch ID of cell where target tree is located
    - `expand::BitVector`: Flag indicating that a patch should be checked for gap expansion 
    (value = 1) or not (value = 0). See `expand_gap()`
    - `previous_species::Vector{Float64}`: Species ID of the last tree to have occupied every cell 
    in the current model grid.
    - `previous_height::Vector{Float64}`: Final height of the last tree to have occupied every cell 
    in the current model grid.
    - `a_height::Float64`: Height in meters of the tree.
    - `id::Int64`: Agent ID of the selected tree/agent
    """
    function rust_impact(agent,
                         rust_infected::Bool,
                         min_symptomic_age::Int64,
                         symptom_prob::Float64,
                         mortality_prob::Float64,
                         gap_maker::Int64,
                         species_ID::Int64,
                         cell::Int64,
                         expand::BitVector,
                         previous_species::Vector{Float64},
                         previous_height::Vector{Float64},
                         a_height::Float64,
                         id::Int64,
                         )

        if rust_infected == true
            agent.rust_infected_age += 1
        end

        #* Check whether tree has become symptomatic
        if agent.rust_infected_age ≥ min_symptomic_age && 
            agent.rust_symptomatic == false
            if rand() < symptom_prob
                agent.rust_symptomatic::Bool = true
            end
        end

        #* If symptomatic chance for tree to die
        if agent.rust_symptomatic == true 
            if rand() < mortality_prob
                if gap_maker[species_ID] == 1
                    expand[cell] = true
                end
    
                ## Record dying tree species and height as a cell list
                previous_species[cell] = species_ID
                previous_height[cell] = a_height
    
                kill_agent!(id, model)
            end    
        end
    end
end