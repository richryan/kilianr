
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kilianr

<!-- badges: start -->
<!-- badges: end -->

The kilianr package provides a system of tools for working with
time-series data, based on the book “Structural Vector Autoregressive
Analysis” by Lutz Kilian and Helmut Lütkepohl. The goal of kilianr is to
reproduce Lutz Kilian’s results from

> Kilian, Lutz. 2009. “Not All Oil Price Shocks Are Alike: Disentangling
> Demand and Supply Shocks in the Crude Oil Market.” American Economic
> Review, 99 (3): 1053–69. DOI: 10.1257/aer.99.3.1053

Kilian’s insights about the identification of oil shocks—and associated
tools—are then used to investigate how oil shocks affect Kern County,
California. In other words, the goal is not creating an R package, per
se; rather, the goal is using Kilian’s tools. This goal has dictated how
the package is constructed. More packages are imported and suggested
than what are strictly necessary to reproduce Kilian’s results.

As a waypoint, I am attempting to publish a description of the package
in [Journal of Open Source Software](https://joss.theoj.org/). The draft
is included in `joss/paper.md`. To help cite work in the paper, I
installed the RStudio Addin [citr](https://github.com/crsh/citr), which
makes accessing entries in a .bib file easier. And the RStudio Addin
[wordcountaddin](https://github.com/benmarwick/wordcountaddin), which
provides a function that counts the number of words in a highlighted
region.

\*\*NOTE: This package uses the package `seasonal`, which is large.

## Philosophy of

## Notes on formatting documentation using `roxygen2`

For formatting documentation, I used
<https://roxygen2.r-lib.org/articles/formatting.html>

The link includes how to make a table.

<!-- ## Installation -->
<!-- You can install the development version of kilianr from [GitHub](https://github.com/) with: -->
<!-- ``` r -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("richryan/kilianr") -->
<!-- ``` -->
<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- # ```{r example} -->
<!-- # library(kilianr) -->
<!-- ## basic example code -->
<!-- # ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
