# Summary statistics for the different variables
final_imputed_data %>% group_by(Waterbody) %>% summarise(sample_count=length(Sample.ID))

final_imputed_data %>% group_by(Trophic_Level) %>% summarise(sample_count=length(Sample.ID))

final_imputed_data %>% group_by(Class) %>% summarise(samp=length(Sample.ID),species = length(unique(Revised_Species)))

