# mIFsim

R package to simulate multiplex immunofluorescence (mIF) data to benchmark spatial statistical methods

# Installing mIFsim to RStudio

To install `mIFsim`, it is required to have `devtools` or `remotes` installed for their `install_github()` function:

```{if (!require("devtools", quietly = TRUE))}
  install.packages("devtools")

devtools::install_github("FridleyLab/mIFsim")
```
