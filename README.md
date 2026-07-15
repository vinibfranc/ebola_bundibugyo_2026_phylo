# ebola_bundibugyo_2026_phylo
Initial treedater dating analysis of EBOV Bundibugyo Ebola outbreak

To generate the Rmarkdown report for the dating analysis:
```
rmarkdown::render(
    input = "scripts/02b_dating.Rmd",
    output_file = "BDBV_2026_dating.html",
    output_dir = "."
)
```