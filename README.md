# Pilgrim_River_Regime

A regime transition model based on an analysis Curry Cunningham prepared for Bristol Bay sockeye. This application is data limited and was done to see if it would fit a 17 year dataset when used with informative priors. The ADF&G network location for this repository is S:\RTS\Reimer\Pilgrim\_sockeye. This analysis was done well into the escapement goal review process and only referred to in the 2024 AYK escapement goal report. It's my hope that it will be updated and considered more seriously during the next cycle. It produced similar estimates of Smsy as the standard analysis but but with better interpretability and more realistic probabilities of acheiving x % or MSY.

-   data: Pilgrim...csv is the raw data received from AYK staff. brood_table.xlsx is a brood table I created to read into jags.

-   quarto: quarto and html files shared with AYK staff. Focused mostly on comparing the results to the simple analysis produced by the escapement goal shiny app.

-   scripts: script I used to develop the model. Not presented staff but there are some interesting results about escapement goal recommendations based on simulated yields or an OYP approach with and without accounting for process error. Would like to explore this correlation more broadly in the future.
