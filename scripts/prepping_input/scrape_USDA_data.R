#code copied from http://bricol.net/downloads/data/PLANTSdatabase/09-02-01PLANTSscrape.R
#Green, 2009

vars <- c('Duration',
          'Growth Habit',
          'Native Status',
          'Federal T/E Status',
          'National Wetland Indicator',
          'Active Growth Period',
          'After Harvest Regrowth Rate',
          'Bloat',
          'C:N Ratio',
          'Coppice Potential',
          'Fall Conspicuous',
          'Fire Resistant',
          'Flower Color',
          'Flower Conspicuous',
          'Foliage Color',
          'Foliage Porosity Summer',
          'Foliage Porosity Winter',
          'Foliage Texture',
          'Fruit/Seed Color',
          'Fruit/Seed Conspicuous',
          'Growth Form',
          'Growth Rate',
          'Height at 20 Years',
          'Height, Mature',
          'Known Allelopath',
          'Leaf Retention',
          'Lifespan',
          'Low Growing Grass',
          'Nitrogen Fixation',
          'Resprout Ability',
          'Shape and Orientation',
          'Toxicity',
          'Adapted to Coarse Textured Soils',
          'Adapted to Fine Textured Soils',
          'Adapted to Medium Textured Soils',
          'Anaerobic Tolerance',
          'CaCO3 Tolerance',
          'Cold Stratification Required',
          'Drought Tolerance',
          'Fertility Requirement',
          'Fire Tolerance',
          'Frost Free Days, Minimum',
          'Hedge Tolerance',
          'Moisture Use',
          'pH, Minimum',
          'pH, Maximum',
          'Planting Density per Acre, Minimum',
          'Planting Density per Acre, Maximum',
          'Precipitation, Minimum',
          'Precipitation, Maximum',
          'Root Depth, Minimum',
          'Salinity Tolerance',
          'Shade Tolerance',
          'Temperature, Minimum',
          'Bloom Period',
          'Commercial Availability',
          'Fruit/Seed Abundance',
          'Fruit/Seed Period Begin',
          'Fruit/Seed Period End',
          'Fruit/Seed Persistence',
          'Propagated by Bare Root',
          'Propagated by Bulb',
          'Propagated by Container',
          'Propagated by Corm',
          'Propagated by Cuttings',
          'Propagated by Seed',
          'Propagated by Sod',
          'Propagated by Sprigs',
          'Propagated by Tubers',
          'Seed per Pound',
          'Seed Spread Rate',
          'Seedling Vigor',
          'Small Grain',
          'Vegetative Spread Rate',
          'Berry/Nut/Seed Product',
          'Christmas Tree Product',
          'Fodder Product',
          'Fuelwood Product',
          'Lumber Product',
          'Naval Store Product',
          'Nursery Stock Product',
          'Palatable Browse Animal',
          'Palatable Graze Animal',
          'Palatable Human',
          'Post Product',
          'Protein Potential',
          'Pulpwood Product',
          'Veneer Product')

index_list <- read.table('./data/raw/09-02-02PLANTSdata.csv', header = TRUE, sep = ' ', colClasses = 'character')

symbols <- index_list$Symbol
data <- matrix('', nrow = length(symbols), ncol = length(vars))

#not sure what this is supposed to do (ended up with NAs)
#for(i in 1:4){
for(i in seq_along(symbols)){
  #  url <- paste('http://plants.usda.gov/java/charProfile?symbol=',
  #               symbols[i], '&format=print', sep = '')
  url <- paste('http://plants.usda.gov/java/charProfile?symbol=',
               symbols[i], sep = '')
  raw <- try(scan(url, what = '', sep = '\n'))
  for(j in seq_along(vars)){
    this_var <- raw[grep(vars[j], raw, perl = TRUE,
                         useBytes = TRUE)[1] + 1]
    this_var <- gsub("<[^<>]*>", "", this_var,
                     perl=TRUE, useBytes = TRUE)
    this_var <- gsub('^ *', '', this_var,
                     perl=TRUE, useBytes = TRUE)
    data[i,j] <- this_var
  } #end loop through variables
} #end loop through species
colnames(data) <- vars
data.out <- cbind(index_list, data)
#write.table(data.out, './data/tidy/USDA_data.out.txt')