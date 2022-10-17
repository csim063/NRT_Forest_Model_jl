# Nga Rakau Taketake Forest Dynamics Model

![Cover_Vid](C:\Users\simpk\Dropbox\Professional\Academic\Uni Auckland\NRT\NRT_Forest_Model_jl\Outputs\Videos\Cover_Vid.gif)

This project implements a spatially explicit individual-based model representing the dynamics of canopy and understory in a New Zealand based forest fragment. The model has been parameterized on data collected primarily in northern New Zealand. See [Morales and Perry (2017)](https://www.sciencedirect.com/science/article/abs/pii/S0304380016306068) and [Brock et *al.* (2020)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.13305) for more details. This version of the model has been implemented in [Julia](https://julialang.org/) making extensive use of the [Agents.jl](https://juliadynamics.github.io/Agents.jl/stable/) package.

*This model has been developed and tested on Windows*.

# Getting started

For information on how to install Julia please see [Download Julia](https://julialang.org/downloads/). We highly recommend using Visual Studio Code with the Julia extension for running and modifying this model.

To install the dependencies used by the model open the Julia REPL and run:

```julia
] add Pkg #To install PKg only if needed remember to exit back to REPL before running next line
using Pkg
```

After running the above code run the `Add_dependencies.jl` script.

## Usage

The model is divided largely into two parts: 

1. Modules which define the functioning of the model and are contained in the `ABM/Modules/` folder.
   - `Initialise.jl`: Module primarily contains the function to define the agents and setup the model for the ABM.
   - `Step.jl`: Module contains the agent and model stepping functions which are the key functions which actually run the model.
   - `Demographic_functions.jl`: Module contains the core functions which enable agents to grow, reproduce, expand their ranges, and die.
   - `Demographic_assignments.jl`: Module contains functions to calculate the demographic parameters (currently just age) to be assigned to trees (agents) during initialization either based on initial height for tree-ferns or DBH for trees.
   - `Disturbance_functions.jl`: Module contains all functions which are used to create disturbances in the model.
   - `Helper_functions.jl`: Script contains two modules both of which contain functions which are largely miscellaneous and do not easily fit into other modules, however are helpful for model scripting. 
2. Scripts to run the model contained in `ABM/`. There are two model running scripts:
   1. `Forest_Model.jl`: This is the primary script of the NRT forest model. This script sets up and runs the model. It can also record and export data from the model run, and visualize the model.
      - You are able to fully define the models behaviors within this script.
   2. `Profile_model.jl`: This scripts runs the NRT forest model at a small scale while running a "sampling" profiler to benchmark the functions called while running the model. For more details regarding how the profiling works see [the Julia profiling documentation](https://docs.julialang.org/en/v1/manual/profile/).

# How to contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Ask a question :interrobang:

Browse the documentation to see if you can find a solution. Still stuck? Open an [issue on GitHub](https://github.com/csim063/NRT_Forest_Model_jl/issues) on GitHub. We'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email](simpkinscraig063@gmail.com).

Please try to include a reproducible example.

### Propose an idea :bulb:

Take a look at the documentation and [issue on GitHub](https://github.com/csim063/NRT_Forest_Model_jl/issues) list to see if it isn't included or suggested yet. If not, please open a new issue!

While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

### Report a bug :bug:

Report it as an [issue on GitHub](https://github.com/csim063/NRT_Forest_Model_jl/issues) so we can fix it. A good bug report makes it easier for us to do so, so please include:

* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug

Again, please try to include a reproducible example.

### Improve the documentation :book:

Good documentation makes all the difference, so your help to improve it is very welcome.

### Pull request process :arrow_up_down:

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [the repo](https://github.com/csim063/NRT_Forest_Model_jl) and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/). Don't forget to pull all new changes before starting to work!

2. Create a new Git branch and use a name that briefly describes the proposed changes.

4. Make your changes:
    * Write your code.
    * Test your code.
    * Document your code (see function documentation above).
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request) and notify a reviewer.
