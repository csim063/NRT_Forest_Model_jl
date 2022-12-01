# Lab book and rough notes

## Substantive changes from NetLogo model

- `capture_gap()`: NetLogo uses previous height as agent is dead at time of calling (NetLogo always has a height of 0). Aditionally, a cone angle of 0 is used in NetLogo meaning only current patch is selected with no further expansion. However, even if the desired 45 degree cone was used this would still differ from Julia implementation which only selects an axis line along which to fall.
  - Further for growth form 1 we use the equation as stated in Morales and Perry 2017, which uses both b2 and b3 jobawa (NetLogo only uses b2). For growth form 2 we assign a DBH. unlike NetLogo.

## Possible areas of improvement

- Currently for in patch dispersal two Poisson distributions are used to represent the split of seeds dispersed under the crown and long distance. This would make more sense and likely be more performant to have as only a single distribution correctly allocated.
- The current implementation the assignment step of `ldd_within()` is very messy and slow. It should be refactored
- In `ldd_within()` the fore loop to assign seeds needs to calculate set differences between all cells in distance and all cells in distance - 1 to get only the cells at distance. It has been suggested that an update to the Agents.jl source code is needed to create a function to just get cells at distance
- A number of functions currently pull from globals rather than having local variables created. This has been shown to be less performant and is likely an easy area for performance gain. Current problem functions are:
  - `herbivore_effect`
  - `thin_regenbank()`
  - `macro_litterfall()`
  - `expand_gap()`
  - `assign_demographics()`
  - Model stepping function
- Currently only the last previous species ID and height is recorded in each cell. This could easily be adapted to record all previous trees if desired.
- Create a kill function to avoid duplicate code used every time a tree is killed
- Currently in `regenerate_patch_bank()` new_saplings is used as the amount of individuals to add to both saplings and seedlings

## Additions

- Added phytothera spread and impacts. The disease model is implemented as an SEI (Susceptible, Exposed, Infected) model. Trees can be defined as susceptible/target species in the demography input file
  Infection spreads (i.e. trees are exposed) in two ways, a global probability of a target tree getting exposed and a local probability that a target tree will get exposed by an infected neighbour in a predefined infectious radius. The more infected neighbour trees the higher the probability (or rather the more tests are done against that probability) a the target tree gets infected. Neighbours may only infect trees if they have been infected longer than a user defined time period to mimic the time taken for trees to become infectious after exposure. Infected trees may become symptomatic after a determined minimum period of time after infection and symptomatic trees may die with a given probability.

- Added grass as an optional factor. If grass is enabled there is a defined probablilty that grass may invade a patch after a tree dies and a gap establishes. If a patch has been invaded by grass any seeds
  landing there only have a probability (rather than a certainty) of establishing into a seedling during dispersal phases. Each dispersal phase and thus each seed entering patch in a step is tested independently.

- Removed previous species for each patch as it is not used and can be accessed directly by data export if needed for analysis.

## Questions

- How realistic is it to empty the regeneration bank (seedlings and saplings) when a tree grows in a gap?
- Should we implement a species specific phytothera spread rate rather than constants.
- Need to assess the radius of surrounding trees assumed to be able to infect trees with phytothera, current it is radius 2 as literature mentioned short range movement

## Possible to do items

- Add impact of edge effect to disease spread. Should be easy enough to just add an extra chance of global infection to edge trees of target species
- Add additional disease impacts resulting from becoming symptomatic.
- Could possibly improve performance by making different agents if there is or isn't disease present as if there is no disease agents require far less parameters
