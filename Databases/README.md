## Problems with downlaoded files

There is sth. weird going on when download the files via curl

``` default
Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  (
README.rmd#14): line 103 did not have 21 elements
```

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

``` r
setwd("/home/td/sshfs_work/clean/MeMoMe")
modelseed_metabolies = read.table("Databases/ModelSEED/compounds.tsv", sep = '\t', head = TRUE)
summary(modelseed_metabolies)
```

    ##       id            abbreviation           name             formula         
    ##  Length:33992       Length:33992       Length:33992       Length:33992      
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##      mass              source            inchikey             charge        
    ##  Length:33992       Length:33992       Length:33992       Min.   :   -2400  
    ##  Class :character   Class :character   Class :character   1st Qu.:      -1  
    ##  Mode  :character   Mode  :character   Mode  :character   Median :       0  
    ##                                                           Mean   :   77959  
    ##                                                           3rd Qu.:       0  
    ##                                                           Max.   :10000000  
    ##     is_core        is_obsolete    linked_compound     is_cofactor
    ##  Min.   :0.0000   Min.   :0.000   Length:33992       Min.   :0   
    ##  1st Qu.:1.0000   1st Qu.:0.000   Class :character   1st Qu.:0   
    ##  Median :1.0000   Median :0.000   Mode  :character   Median :0   
    ##  Mean   :0.8147   Mean   :0.001                      Mean   :0   
    ##  3rd Qu.:1.0000   3rd Qu.:0.000                      3rd Qu.:0   
    ##  Max.   :1.0000   Max.   :1.000                      Max.   :0   
    ##      deltag           deltagerr            pka                pkb           
    ##  Min.   :  -11540   Min.   :       0   Length:33992       Length:33992      
    ##  1st Qu.:     -66   1st Qu.:       2   Class :character   Class :character  
    ##  Median :     159   Median :       6   Mode  :character   Mode  :character  
    ##  Mean   : 4172125   Mean   : 4172160                                        
    ##  3rd Qu.:10000000   3rd Qu.:10000000                                        
    ##  Max.   :10000000   Max.   :10000000                                        
    ##  abstract_compound  comprised_of         aliases             smiles         
    ##  Length:33992       Length:33992       Length:33992       Length:33992      
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##     notes          
    ##  Length:33992      
    ##  Class :character  
    ##  Mode  :character  
    ##                    
    ##                    
    ## 
