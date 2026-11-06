# agentcrabr

Agent-based model (ABM) for crab population dynamics, implemented in R using the
ABM package, with an interactive Shiny dashboard for visualization.

## Features

- Crab agents with states: foraging, resting, molting, dead
- Energy dynamics: gain during foraging, cost during resting/molting
- Growth via molting
- Mortality from starvation
- Interactive Shiny dashboard with population dynamics, state distribution, and more

## Installation

`
# Install dependencies
install.packages(c("ABM", "shiny", "ggplot2", "dplyr", "tidyr"))
`

## Usage

`
# Run the Shiny app
shiny::runApp("app.R")
`

Or from the command line:
`ash
Rscript -e "shiny::runApp('app.R')"
`

## Project Structure

- R/model.R - Crab ABM model definition
- pp.R - Shiny dashboard application
- inst/ - Additional resources
- www/ - Web assets

## Model Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| n_crabs | 50 | Initial number of crabs |
| energy_init | 10 | Initial energy level |
| energy_max | 20 | Maximum energy capacity |
| energy_gain_foraging | 0.3 | Energy gained per foraging event |
| energy_cost_resting | 0.05 | Energy cost per resting period |
| energy_cost_molting | 3 | Energy cost of molting |
| molt_growth | 0.1 | Size increase per molt |
| foraging_rate | 0.5 | Rate of transitioning to foraging |
| resting_rate | 0.3 | Rate of transitioning to resting |
| molting_rate | 0.05 | Rate of transitioning to molting |
| mortality_rate | 0.02 | Base mortality rate |

## References

- ABM R Package: https://cran.r-project.org/package=ABM
- Shiny: https://shiny.posit.co/
