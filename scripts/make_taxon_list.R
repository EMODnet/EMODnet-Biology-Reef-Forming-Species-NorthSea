#install the required packages (comment if needed)
# install.packages("jsonlite", repos="http://cran.r-project.org")
# install.packages("httr")

#Use the libraries
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)
library(dplyr)
library(tidyr)
library(purrr)

#Fill in the AphiaID you need
AphiaID <- c(130867, 131495 , 140467, 140658)

#Build the URL to get the data from
url <- lapply(AphiaID, function(x) sprintf("https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/%d", x))

#Get the actual data from the URL
classificationTree <- lapply(url, function(x) (fromJSON(x)))
# names(classificationTree) <- AphiaID

#Walk the classification tree

ranks <- list()

for(j in 1:length(AphiaID)){
  
  currentTreeItem = classificationTree[[j]]
  currentAphiaID = AphiaID[j]
  species_df = tibble()
  i = 1
  
  while (!is.null(currentTreeItem )) {
    
    species_df_temp = 
      tibble(
        AphiaID_species = currentAphiaID,
        aphiaid = (currentTreeItem$AphiaID),
        scientificname = print(currentTreeItem$scientificname),
        rank = (currentTreeItem$rank)
      )
    
    if(i == 1) {
      species_df <- species_df_temp
    }  else {
      species_df <- species_df %>% bind_rows(species_df_temp)
    }
    
    i = i + 1
    #Get next item in the tree
    currentTreeItem <- currentTreeItem$child;
  }

  ranks[[j]] <- species_df
  
    
}

taxon <- lapply(ranks,
       function(x){
         pivot_wider(x, id_cols = AphiaID_species, names_from = rank, values_from = scientificname)
       }
) %>% 
  bind_rows() %>%
  dplyr::select(
    AphiaID_species,
    Phylum,
    Class,
    Genus,
    Species
  )

rm(ranks, species_df_temp, species_df, url, classificationTree)


