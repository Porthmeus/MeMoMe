``` r
setwd("/home/td/sshfs_work/clean/MeMoMe")
bigg_metabolies = read.table("Databases/BiGG/bigg_models_metabolites.tsv", sep = '\t', head = TRUE)
summary(bigg_metabolies)
```

    ##    bigg_id          universal_bigg_id      name            model_list       
    ##  Length:15724       Length:15724       Length:15724       Length:15724      
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##  database_links     old_bigg_ids      
    ##  Length:15724       Length:15724      
    ##  Class :character   Class :character  
    ##  Mode  :character   Mode  :character
