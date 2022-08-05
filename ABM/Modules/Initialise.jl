"""
Module contains all of the functions required to setup the model at the start
of a run
"""
module Setup
using Agents
using Random
    #% Define agents
    #! !!!! BEST APPROACH IS PROBABLY TO HAVE ADULT TREES AS AGENTS AND 
    #! EVERYTHING ELSE INCLUDING SAPLINGS/SEEDLINGS AS PATCHE PROPERTIES
    #! Not that i have placed defining agents first as agents are required to
    #! actually initialise a model instance
    ## TODO: Currently just dummy agents to try get model working
    # TODO: probably will need to define different agents for different things
    #? Maybe a dict of species key value pairs could be useful
    @agent Trees GridAgent{2} begin 
        species_ID::Int
        height::Int
        dbh::Int
        age::Int
    end

    #% Define world
    function forest_model(;
        griddims = (50, 50),
        site_df = site_df,
        seed = 999,)
        
        #? Could we define space as many dimensional and just add the properties like height
        #? and species present to this
        space = GridSpaceSingle(griddims; periodic = false);
        rng = MersenneTwister(seed)

        model = ABM(Trees, space; 
            rng,
            scheduler = Schedulers.Randomly())

        ## Populate the world with adult tree agents
        grid = collect(positions(model))
        num_positions = prod(griddims)

        #Make for loop that samples a proportion of space and allocates each species
        #? Maybe find a way to get one tree per cell and then just have species set to a list with proportions set by proportion
        used_positions = []::Matrix{Tuple{Int64, Int64}}
        for species in 1:nrow(site_df)
            #Create sample of
        end
        
        return model
    end

    #% Define patches

end