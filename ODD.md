# Nga Rakau Taketake Forest Dynamics Model Description

*Note that this model description is heavily derived from the original ODD document written for [Morales and Perry (2017)](https://www.sciencedirect.com/science/article/abs/pii/S0304380016306068) and has only been updated to reflect the changes made to the model. As such it should not be used as is in any publication.*

## Overview

### Model description

The model description below follows the Overview, Design concepts, and Details (ODD) protocol described by Grimm et al. (2006, 2010); this approach is designed to support comprehensive and transparent model description and communication. To ensure that this model description is complete some information is repeated from the general description presented in the main document.

### Purpose

The main purpose of our model was to enable us to evaluate and/or assess the long-term viability of kauri or pohutakawa forest fragments (though the model is generalised for any NZ forest) under different management scenarios (including no management), disease spreads, and climatic conditions.

### Entities, state variables, and scales

The basic unit represented in the model was an individual tree, with individuals divided into three distinct life-stages: seedlings (individuals < approximately 50 cm in height), saplings (individuals > approximately 50 cm in height) and adults (Figure 1) ([Morales and Perry (2017)](https://www.sciencedirect.com/science/article/abs/pii/S0304380016306068) section 2.1; model description). Each individual adult tree was characterized by the following state variables: spatial location (on a discrete-space grid or lattice), species identity, diameter at breast height (dbh), height, age, growth form (either tree or tree fern), and disease infection status (for both phytothera and rust). Individual seedlings and saplings did not have an explicit spatial location; instead these life-stages were represented by their abundance in a grid cell.  Each species was characterized by the following parameters: shade tolerance, reproduction age, regeneration height, seed production, dispersal, gap maker, suppression tolerance, suppression mortality rate, palatability (herbivory), edge response, neighbourhood dispersal and in-fragment dispersal (each explained further below). Seedlings and saplings were characterized by their survival, and seedlings by the rate at which they transition to the sapling class.

The model was initialized and parameterized using the same data as used in [Morales and Perry (2017)](https://www.sciencedirect.com/science/article/abs/pii/S0304380016306068), that being data collected from unfragmented forest (continuous fenced forest) and from the published literature. The study area was located near Cambridge in the Waikato region, North Island, New Zealand. The unfragmented forests were located at Te Miro Scenic Reserve, Maungatautari Ecological Island North and Maungatautari Ecological Island South (Morales et al., 2016). The Maungatautari Ecological Island is surrounded by a predator-proof fence (Burns et al., 2012) and Te Miro is subject to active control of exotic mammalian pests (Morales et al., 2016).

The model was initialized with a fixed number of seedlings, saplings and trees specific for each species, with adult abundances differing for each site category (unfragmented forest, fenced, unfenced fragments). Each adult tree was allocated a diameter at breast height (dbh) from which its age and height were estimated based on standard allometric curves. The initial distribution of dbh values was simulated from a log-normal distribution conditioned on size-frequency information from each site. Thus, each time the model is run the initial size-frequency structure of the forest will differ, but the composition and structure of the forest will be similar.

One time-step or ‘tick’ (the model’s fundamental temporal unit) represents one year and, unless otherwise noted. The model is grid-based and only one adult tree can occupy any one grid cell, with the grid cells representing a 4 × 4 m cell by default, though this can be adjusted by the user. We choose this grid resolution because it is the approximate canopy size of an individual B. tawa tree and is close to the mean inter-tree distance recorded in the field. The model is spatially explicit and can simulate the dynamics of each individual adult tree in forest areas up to a size of ca. 16 ha; larger areas can be simulated but most of the forest fragments in the study area are very small (Waikato Regional Council, 2009).

### Process overview and scheduling
In each tick the model two types of steps are run, initially processes needing to be run for all agents, followed by processes run globally across the whole model, these processes are a series of ecological routines run sequentially for each type of step. The processes run for all agents are initial demography, growth, within patch seed dispersal, herbivory, disease, and background mortality. Global processes are gap formation via disturbance (such as fire), beyond patch seed dispersal, gap expansion and reduction, restoration, regeneration (recruitment), seedling/sapling demography and ENSO severity calculation.

Before running the model, the user can define different parameterisations and management scenarios. The model can be run with or without management, with or without disease, with different levels of disturbance or herbivory, with or without grass, and with or without impacts from weather altered by ENSO cycles. The user can also define the number of ticks (years) to run the model for.

At initialization the species abundances, and the dbh and age of each individual adult tree depend on the scenario being assessed (unfragmented forest, fenced and unfenced fragment). The model starts with a forest with the same number of seedlings (10 seedlings [individuals < approximately 50 cm in height] per species per grid cell) regardless of the scenario (up to 50 cm in height) and saplings (one sapling [individuals > approximately 50 cm in height) per species per grid cell). The number of seedlings was highly variable across the different sites from where we had collected parameterisation information. We choose to use the average of the total number of seedlings, which corresponded to approximately 8306 individuals per ha or 8.3 individuals per cell; we rounded this up to 10 individuals per grid cell. 

The transition from seedling to adult tree has two stages. The first stage corresponds to seedlings that form the seedling bank.  A species-specific proportion of seedlings transition to saplings each year.  Each year individuals in the sapling stage can transition (via a lottery process weighted by the light levels in each grid cell) to adult trees if there is not already an adult tree present in the grid cell. In both life-history stages, if the individuals survive but do not transition to the next stage they remain as part of the seedling or sapling bank until either they successfully move to the next stage or they die. A species-specific proportion of seedlings and saplings dies each year.  

Natural regeneration, as simulated, starts with dispersal processes occurring at three distinct scales: neighbourhood dispersal (‘under parent tree’ dispersal), in-fragment dispersal (movement of seeds away from parent individuals but within the fragment), and long-distance dispersal from beyond the simulated fragment.  In order to reduce the computational burden, and as is common in forest models, the fate of individual seeds is not modelled, and instead established seedlings are ‘dispersed’. Both local and long-distance dispersal processes depend on the ability of each species to disperse by vectors such as wind and birds. If restoration activities are simulated, the planting of saplings is represented by adding a user-determined number of saplings of each species to the sapling bank.

Mortality consists of three broad components: an annual background mortality rate, suppression arising from slow growth or competition (due to neighbouring individuals competing for gradients such as light), and disease. The annual background mortality rate is calculated following Shugart (1984; details below) and is age-independent and impacted by ENSO severity (if this option is chosen), with individuals having an approximately 2% chance of surviving to their maximum possible age.  If an individual’s growth (averaged over the previous five years) is less than the species critical growth (a proportion of the optimal growth for an individual at a given age, with lower values indicating increasing shade-tolerance or ability to persist through long periods of suppression), it is assumed that the individual is suppressed and hence suffering elevated stress and is more vulnerable to disease and so forth (what Keane et al. 2001 term ‘growth-dependent mortality’). Suppressed individuals die at a species-specific rate, which is higher than background mortality. Diseases are modelled using a SEI (Susceptible, Exposed, Infected). There are two types of disease modelled, phytophera (soil-borne) and rust (air-borne). Trees are defined as susceptible/target species to each disease type individually. Infection spreads (i.e. trees are exposed) in two ways, a global probability of a target tree getting exposed and a local probability that a target tree will get exposed by an infected neighbour in a predefined infectious radius. The more infected neighbour trees the higher the probability (or rather the more tests are done against that probability) a the target tree gets infected. Neighbours may only infect trees if they have been infected longer than a user defined time period to mimic the time taken for trees to become infectious after exposure. Infected trees may become symptomatic after a determined minimum period of time after infection and symptomatic trees may die with a given probability. If a dead tree is a ‘gap-maker’ species, a gap is formed shaped as a cone of a length proportional to the height of the dead tree in a random direction from the tree. All smaller (shorter) adult individuals in patches falling under that cone are assumed to have been hit and damaged by the falling tree and die.

Each canopy-dominant tree grows following the equations described by Botkin et al. (1972) and that underpin the JABOWA gap model.  Botkin’s approach calculates the optimal diameter growth rate for an individual of a given species and size and then down-weights it based on competitive interactions and environmental conditions. An individual tree’s height is calculated following the allometric approach described by Botkin et al. (1972) and Botkin (1993); full details are provided below. In our model the diameter increment is reduced by competition and an edge effect, but; environmental conditions such as soils are assumed to be uniform. As we were simulating small forest fragments we assumed that soil properties did not vary within such small areas. A competitive effect is imposed if a tree’s height is less than any of the neighbouring trees’ heights. If edge effects are represented then a penalty is applied to the growth increment depending on the individual’s location in the plot and how well the species can tolerate the edge environment. The edge effect is represented by a species-specific parameter based on published information describing how each species included in the model responds to edge environments.

### Design concepts

#### Basic principles

The guiding objective of this model was to simulate the long-term dynamics of native kauri and pohutakawa forest fragments in northern New Zealand. In order to simulate different management activities such as fencing and herbivore control, several Boolean switches were added: herbivory, edge effects, the active planting of seedlings (restoration activity), diseases, and weather impacts. The model represents fundamental demographic processes such as seed dispersal, reproduction, and mortality. Spatio-temporal variability in soils were not explicitly included, however a generic gradient impacting tree growth has been implemented and may be parameterised to match a soil nutrient gradient. Individual species demographic rates are parameterized using data collected specifically for this purpose, supplemented by empirical data from other studies; in particular, maximum height, shade tolerance and maximum age were determined from the available literature. 

#### Emergenence

Forest dynamics and population structures emerge from the interaction of individuals across the three life-stages that we considered (seedling, sapling and adult tree).

#### Sensing

Individuals “know” their life cycle stage (seedling, sapling, and adult tree) and what species they belong to; they also “know” their age, height and their diameter.

#### Interaction

Competition for light is the fundamental interaction represented in the model. When the amount of light that a given individual receives is insufficient to sustain growth (a species-specific value) then the individual is considered to be suppressed. The point at which this suppression occurs depends on each species shade-tolerance. The period over which a suppressed individual can remain in this state without experiencing elevated mortality is species specific. A second generic gradient is also included in the model, which impacts growth and is assumed to represent a soil nutrient gradient. This gradient is based on values in each cell drawn from a random uniform distribution between 0 and 1. This value is used to in association with an exponential distribution to calculate the probability density for that gradient value on that distribution. This probability density is scaled between 0 and 1 and then subtracted from 1 to ensure cells with high values (i.e. higher resources) have higher competitive values and thus higher growth. The generic gradient is used to assess what species grow in a gap, where it is used as the mean of a truncated normal distribution used to draw a weighting value for each species at each step. This weighting value is combined with the shade weighting value with the mean used as the overall weight to be used when selecting the species in the lottery function.

#### Stochasticity

Nearly all model components are stochastic.

### Details

#### Initialisation

Species are initially randomly assigned to each of the grid cells. The initial number of seedlings and saplings present in each cell was defined using empirical data collected from each of the study sites. Adult trees’ dbh are generated randomly from a log-normal distribution with mean and standard deviation estimated from field data. The data used in this initialization process come from three unfragmented forest sites. Tree heights and age are calculated from the initial diameter (see submodel section).  By consequence, a different forest is generated each time the model is set up, although its overall statistical properties are the same. The remainder of the parameters that comprise the species demography, such as shade-tolerance, maximum height, among others, are fixed across a species lifetime and are constant over all simulations. The pervasiveness of the edge effect is controlled by the parameter e2, which was always 0.3.

#### Input data

The model uses input data to represent the probability of transition between ENSO states. This data was derived using the 100 years of NIWA data taken from https://niwa.co.nz/climate/information-and-resources/elnino.

### Submodels

#### Initial demographics

Patches are initialized as described above using empirical data collected for this purpose.

Each individual’s growth is represented by an annual diameter increment ($\Delta D$; m/yr), which is calculated using the formula provided by Botkin et al. (1972) for the JABOWA gap model. The gap model approach has not been used frequently in New Zealand (but see DeVelice, 1988 and Hall and Hollinger, 2000). However, the standard approach yields an optimal growth curve for a given species.

$$\Delta D = \frac{(dbh \times g \times (1-(dbh \times H)/(D_{max} \times H_{max})))}{(2.74 + 3 \times b_2 \times dbh - 4 \times b_3 \times dbh^2)}$$

where: dbh is diameter at breast height (m), H is tree height in meters, Hmax is maximum attainable height (m) and Dmax is maximum attainable diameter (dbh), which are species-specific. $b_2$ and $b_3$ are species-specific allometric constants. Two-thirds of the maximum diameter ($D_{max}$) is reached at 50 % of the maximum age. $b_2$ and $b_3$ are calculated using the formulae given by Botkin et al. (1972):

$$b_2 = 2(H_{max} - 1.37)/D_{max}$$
$$b_3 = (H_{max} - 1.37)/D_{max}^2$$

where: 1.37 is the height (m) of the smallest adult trees represented in the model, $H_{max}$ is maximum height and $D_{max}$ is the maximum diameter (dbh) a given species can attain.

Diameter increment ($\Delta D_c$) is weighted by a competitive penalty (value from 0 to 1) and edge effect (value from 0 to 1) as follows (assuming edge and neighbour effects are multiplicative):

$$\Delta D_c = D \times CP \times EP$$

where: $\Delta D_c$ is the current diameter increment after the correction, CP is the competitive penalty and EP is edge penalty.

#### Height

The height (in meters) of a tree of known dbh is calculated using the formula provided by Kerr and Smith (1955), which is the standard approach for gap models and their derivatives:

$$H = 1.37 + (b_2 \times dbh) + (b_3 \times dbh^2)$$

where 1.37 is the height in meters of saplings when they transition to adults.

#### Restoration planting

Restoration planting is simulated by adding saplings to the sapling bank at a species-specific amount at a specific frequency (e.g. yearly, every five years, etc.).

#### Dispersal

Each year in each grid cell a number of seedlings are produced, following a Poisson deviate specified by a seed production parameter (seed-prod) (seed production is set to the same value for all the species). Dispersal is comprised of three separate processes: neighbourhood dispersal, in-fragment dispersal and long-distance dispersal. 

##### 1. Neighbourhood dispersal

The crown size ($Cr_{size}$) of each adult tree is computed using their dbh to determine which cells are under the trees following the simplified allometric relationship from SORTIE-NZ (Kunstler et al., 2013; 2009):

$$Cr_{size} = 0.284 \times (d^{0.654})$$

$Cr_{size}$ corresponds to the estimation of crown size of an adult tree based on its dbh in cm (d).

Then a fraction of those seedlings is distributed across the neighbouring cells of the parent trees (the proportion of seedlings is controlled by ldd-dispersal-frac)

##### 2. In-fragment dispersal

A fraction of seedlings (ldd-dispersal-frac) are dispersed within the fragment but beyond their parents. Each such seedling is dispersed a negative exponential distance at a random direction, with a species-specific mean dispersal distance (ldd-dispersal-dist).  The dispersed individual is then added to the seedling bank in the appropriate grid cell.

##### 3. Long-distance dispersal

Each grid cell has a chance of receiving a seedling from outside the simulated fragment. The probability of receiving a seedling is the larger of the species proportional abundance in the grid or a minimum chance of 5%. The seedling actively dispersed by each species is set by the parameter external-species. This parameter is species-specific because some species are more likely to be dispersed by frugivorous birds than others. Then, the number of seedlings for each species that is dispersed is calculated as a binomial process and those seedlings are allocated randomly across the grid. This approach makes the assumption that the composition of the modelled area mimics that of the broader landscape in which it is embedded. Once the aggregate seed rain is calculated, seedlings are dispersed at random across the simulated landscape by adding a seedling to the seedling bank in the appropriate grid cell.

#### Herbivory

We considered herbivory as a component of mortality as described in the Mortality Section below; we mention it here to preserve the sequential order of process scheduling in the model.

#### Disease

Two types of diseases are modelled, phytophthora and rust. Phytophthora is a soil-borne fungus that causes root rot and is a major cause of mortality in New Zealand forests. Rust is a fungal disease that causes leaf spots and can cause mortality in some species. Both diseases only affect species listed as susceptible. The disease model is implemented as an SEI (Susceptible, Exposed, Infected) model. Infection spreads (i.e. trees are exposed) in two ways, a global probability of a target tree getting exposed (rust only spreads via this method) and a local probability that a target tree will get exposed by an infected neighbour in a predefined infectious radius. The more infected neighbour trees the higher the probability (or rather the more tests are done against that probability) a target tree gets infected. Neighbours may only infect trees if they have been infected longer than a user defined time period to mimic the time taken for trees to become infectious after exposure. Infected trees may become symptomatic after a determined minimum period of time after infection and symptomatic trees may die with a given probability.

#### Mortality

##### Juvenile mortality

Each year seedlings and saplings suffer background mortality at a species-specific rate, with the number of individuals in each life-stage dying drawn from a binomial distribution (Figure A.1.4).  Total seedling mortality is given by the total number of individuals minus the proportion of individuals that survive minus the individuals that make the transition to the next life-stage (sapling). Sapling mortality is represented by the total number of saplings minus the proportion of individuals that survive. Remaining saplings are susceptible to mortality from herbivory at a user defined rate.

##### Adult mortality

Adult mortality is the combined result of background mortality, disease, herbivory and disturbance.

##### Background mortality

Annual tree (adult) mortality follows the formula of Shugart (1984):

$$M_{rate} = \frac{4}{m_{age}}$$

$M_{rate}$ yield a background mortality rate such that 98% of individuals die before reaching the species maximum age ($m_{age}$).

#### Gap making

If the dying individual belongs to a gap-maker species then it creates a canopy gap potentially extending beyond its own grid cell on its death. The newly generated gap has a conic shape, with all individuals in affected grid cells suffering mortality.

#### Regeneration

New seedlings are added to each grid cell’s seedling bank each year as described above; saplings may be added either due to seedlings transitioning or by the restoration planting process. Each sapling has a set annual probability of 25% of capturing a gap (i.e. transitioning to adult) if there is no adult tree in the cell; with lottery competition weighted on basis of local light environment if there are multiple species.

#### Edge effects

The edge effect ($E_{eff}$) is represented by the following function:

$$E_{eff} = e_1 \times \exp(-e_2 \times d) + e_0$$

where: e1 = 1, and e2 controls the distance (d) over which the edge effect extends into the forest fragment.

As e2 increases, the distance over which the edge has an effect reduces and vice versa.  As edge effect is calculated for each patch, distance (d) represents the distance from the edge to a given patch. The edge effect is the same on all edges of the simulated fragment.

##### Edge penalty

The edge penalty ranges from 0 to 1 (i.e. the effect of the edge on a species growth-rate increasing mortality risk) and is calculated by multiplying the edge effect relative to distance by the ability of the species to persist in edge environments:

$$EP = (1 - (E_{eff} \times (1 - E_{resp})))$$

where: edge effect ($E_{eff}$) is the edge effect relative to the distance from the edge, $E_{resp}$ ranges from 0 to 1 and represents how well a species is adapted to edge environments.

#### Competition
The strength (scaled from 0 to 1) of competition that each individual is experiencing is represented phenomenologically and is based on an individual’s height relative to its neighbours, which is a proxy for local light availability [based on Dislich et al., (2009)]:

$$CP = \alpha (\frac{H}{H_{nhb}})^{0.5}$$

where: CP is the competitive penalty, α is a constant, which we set to 1.6), H is height from equation 1 and Hnhb is the height of the surrounding trees. 

This competition from light is combined with competition from an as of now generic gradient, with the mean taken as the impact on growth. The generic gradient based on values in each cell drawn from a random uniform distribution between 0 and 1. This value is used to in association with an exponential distribution to calculate the probability density for that gradient value on that distribution. This probability density is scaled between 0 and 1 and then subtracted from 1 to ensure cells with high values (i.e. higher resources) have higher competitive values remembering lower values result in lower growth

#### Weather

Weather may be selected and is used to alter background mortality. The impact of weather is based on a value drawn from a uniform distribution with a minimum of zero and a user defined maximimum. The extent of the distribution is truncated based on the current ENSO state with the lower half being used for weaker "LNL" and "ENL" states and the upper half used for stronger "LN" and "EN" states, neutral "N" ENSO state has an effect value of 0.0. The value is converted to a negative for La Nina states meaning that La Nina reduces mortality and is thus assumed to benefit trees overall. The value drawn from the distribution is added to baseline mortality, so positive values will result in higher mortality rates and lower values will result in decreased mortaility. ENSO conditions have five states (La Nina "LN", Neutral "N", El Nino "EN", La Nina weak "LNL", El Nino weak "ENL"). The ENSO state is updated at the start of each year and is based transition probabilities. These transition probabilities were calculated using 100 years of data taken from NIWA.

## References

- Botkin, D.B., 1993. Forest dynamics: an ecological model. Oxford University Press.
- Botkin, D.B., Janak, J.F., Wallis, J.R., 1972. Some ecological consequences of a computer model of forest growth. J. Ecol. 60, 849–872.
- Burns, B., Innes, J., Day, T., 2012. The use and potential of pest-proof fencing for ecosystem restoration and fauna conservation in New Zealand, in: Somers, M.J., Hayward, M. (Eds.), Fencing for conservation. Springer New York, pp. 65–90.
- Develice, R.L., 1988. Test of a forest dynamics simulator in New Zealand. N. Z. J. Bot. 26, 387–392. 
- Dislich, C., Günter, S., Homeier, J., Schröder, B., Huth, A., 2009. Simulating forest dynamics of a tropical montane forest in South Ecuador. Erdkunde 63, 347–364.
- Grimm, V., Berger, U., DeAngelis, D.L., Polhill, J.G., Giske, J., Railsback, S.F., 2010. The ODD protocol: A review and first update. Ecol. Model. 221, 2760–2768.
- Hall, G.M.J., Hollinger, D.Y., 2000. Simulating New Zealand forest dynamics with a generalized temperate forest gap model. Ecol. Appl. 10, 115–130.
- Ker, J.W., Smith, J.H.G., 1955. Advantages of the parabolic expression of height-diameter relationships. For. Chron. 31, 236–246.
- Kunstler, G., Coomes, D.A., Canham, C.D., 2009. Size-dependence of growth and mortality influence the shade tolerance of trees in a lowland temperate rain forest. J. Ecol. 97, 685–695. 
- Morales, N.S., Perry, G.L.W., Burns, B.R., 2016. Fencing is not enough to reinstate regeneration: evidence from a large fruited canopy tree Beilschmiedia tawa. For. Ecol. Manag. 376, 36–44.
- Shugart, H.H., 1984. A theory of forest dynamics: the ecological implications of forest succession models. Springer-Verlag.
- Waikato Regional Council, 2009. Environmental information. Forest fragmentation. Key points. [WWW Document]. URL http://www.waikatoregion.govt.nz/Environment/Environmental-information/Environmental-indicators/Biodiversity/forest-fragmentation-report-card/ (accessed 8.18.16).

