#Code to reproduce forrest plots in Ward et al. 2024 https://www.medrxiv.org/content/10.1101/2024.12.11.24318702v1
# Depednancies
library(ggplot2)
library(dplyr)
library(gridExtra)
library(readxl)
library(tidyr)
library(patchwork)
################################################################################################

# assume current working directory is the repository root
#Reproduction number
data <- read_excel(file.path("data", "H5N1pptdat.xlsx"), sheet = "RO ")

# Clean column names to remove unwanted spaces or characters
colnames(data) <- trimws(colnames(data))

custom_order <- c(
  "Biggerstaff et al. 2014 - Seasonal influenza (community settings)",
  "Larrauri Camara et al. 2012 - 2009 Pandemic (Spain)",
  "Lee et al. 2021 - 2009 Pandemic (Korea)",
  "Biggerstaff et al. 2014 - 2009 Pandemic (confined settings)",
  "Biggerstaff et al. 2014 - 2009 Pandemic (community settings)",
  "Biggerstaff et al. 2014 - 1968 Pandemic (confined settings)",
  "Biggerstaff et al. 2014 - 1968 Pandemic (community settings)",
  "Biggerstaff et al. 2014 - 1957 Pandemic (community settings)",
  "Biggerstaff et al. 2014 - 1918 Pandemic (confined settings)",
  "Biggerstaff et al. 2014 - 1918 Pandemic (community settings)",
  "Ferguson et al. 2004 - Epidemic in Asian bird populations 2004",
  "Bettencourt and Ribeiro 2008 - Indonesia 2004–2006",
  "Bettencourt and Ribeiro 2008 - Vietnam 2004–2006",
  "Aditama et al. 2012 - Indonesia 2005–2009",
  "Yang et al 2007 - Indonesia 2006",
  "Saucedo et al. 2019 - Family cluster data",
  "Estimated - 2024 US H5N1 outbreak (Scenario 2)",
  "Estimated - 2024 US H5N1 outbreak (Scenario 1)"
)


# Prepare the data for plotting
data <- data %>%
  mutate(
    Plot_Label = paste(Author, "-", Outbreak),  # Combine Author and Outbreak for Y-axis labels
    Plot_Label = factor(Plot_Label, levels = custom_order),  # Apply custom order
    Reproduction_number = as.numeric(Reproduction_number),  # Ensure numeric format
    `95%_lower` = as.numeric(`95%_lower`),
    `95%_higher` = as.numeric(`95%_higher`),
    IQR_lower = as.numeric(IQR_lower),
    IQR_higher = as.numeric(IQR_higher)
    #`95%_CrI_lower` = as.numeric(`95%_CrI_lower`),
    # `95%_CrI_higher` = as.numeric(`95%_CrI_higher`)
  )

# Define subtype colors
subtype_colors <- c("H1N1" = "darkgreen", "H2N2" = "darkblue", 
                    "H3N2" = "orange", "H1N1pdm09" = "grey", 
                    "H5N1" = "red", "Seasonal" = "magenta")


#Forest plot
R0 <- ggplot(data, aes(x = Reproduction_number, y = Plot_Label)) +
  # Add error bars for 95% confidence intervals, with linetype and colored by subtype
  geom_errorbarh(
    aes(xmin = `95%_lower`, xmax = `95%_higher`, color = Subtype, linetype = "95% CI"),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  # Add error bars for IQR, with linetype and colored by subtype
  geom_errorbarh(
    aes(xmin = IQR_lower, xmax = IQR_higher, color = Subtype, linetype = "IQR"),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  # Add points for R₀ values, colored by subtype
  geom_point(aes(color = Subtype), size = 5, na.rm = TRUE) +
  # Customize colors for subtypes
  scale_color_manual(
    values = subtype_colors,
    name = "Subtype"  # Legend title for subtype colors
  ) +
  # Customize linetypes
  scale_linetype_manual(
    values = c("95% CI" = "longdash", "IQR" = "solid"),
    name = "Interval Type"  # Title for the line type legend
  ) +
  # Customize labels and titles
  labs(
    title = "",
    x = "Reproduction Number",
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),         # General text size
    axis.text.y = element_text(size = 16, face = "bold"), # Y-axis text size
    axis.text.x = element_text(size = 16, face = "bold"), # Bold X-axis text
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(1.5, "cm")      # Adjust the overall size of legend keys
  )

ggsave(
  filename = file.path("plots", "R0_review.png"),
  plot = R0,
  width = 14,
  height = 12,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

################################################################################################
#incubation period

data <- read_excel(file.path("data", "H5N1pptdat.xlsx"),  sheet = "Inc")
colnames(data) <- trimws(colnames(data))

# Prepare data for plotting
custom_order <- c(
  "Lessler et al. 2009 - Influenza B",
  "Lessler et al. 2009 - Influenza A",
  "Cao et al. 2009 - 2009 Pandemic (China)",
  "Shen et al. 2012 - 2009 Pandemic (China)",
  "Wang et al. 2012 - 2009 Pandemic (Luoyang  school)",
  "Lessler et al. 2009 - 2009 Pandemic (New York school)",
  "Tom et al. 2010 - 2009  Pandemic (UK)",
  "Nishiura and Inaba 2011 - 2009 Pandemic (Japan)",
  "Nishiura 2007 - 1918 Pandemic (Australia)",
  "Canini and Carratt 2011 - Volunteer challenge",
  "Carrat et al. 2008 - Volunteer challenge",
  "Ferguson et al. 2005 - outbreak  aboard a commercial airliner",
  "Huai et al. 2008 - China 1997-2008",
  "Beigel et al. 2005 - Thailand and Vietnam in 2004",
  "Yang et al. 2007 - Indonesia 2006",
  "Oner et al. 2006 - Turkey 2006",
  "Cowling et al. 2013 - China 2013"
)
# Prepare the data for plotting, removing NA in estimate type
data <- data %>%
  mutate(
    Plot_Label = paste(Author, "-", `Outbreak (Subtype)`),  # Combine Author and Outbreak for Y-axis labels
    Plot_Label = factor(Plot_Label, levels = custom_order),  # Apply custom order
    Incubation_period = as.numeric(Incubation_period),  # Ensure numeric format
    `95%_lower` = as.numeric(`95%_lower`),
    `95%_higher` = as.numeric(`95%_higher`),
    IQR_lower = as.numeric(IQR_lower),
    IQR_higher = as.numeric(IQR_higher),
    range_lower = as.numeric(range_lower),
    range_higher = as.numeric(range_higher)
  )

interval_linetypes <- c(
  "95% CI" = "longdash",  # Dashed line for 95% CI
  "IQR" = "solid",      # Solid line for IQR
  "Range" = "solid"    # Dotted line for Range
)

subtype_colors <- c(
  "H1N1" = "darkgreen",
  "H2N2" = "darkblue",
  "H3N2" = "orange",
  "H1N1pdm09" = "grey",
  "H5N1" = "red",
  "Seasonal" = "magenta",
  "Influenza A" = "purple",
  "Influenza B" = "brown"
)

estimate_shapes <- c("Mean" = 17, "Median" = 16)  # Triangle for Mean, Circle for Median

# Create the forest plot for incubation periods
inc <- ggplot(data, aes(x = Incubation_period, y = Plot_Label)) +
  # Add error bars for 95% confidence intervals
  geom_errorbarh(
    aes(xmin = `95%_lower`, xmax = `95%_higher`, linetype = "95% CI", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  # Add error bars for IQR
  geom_errorbarh(
    aes(xmin = IQR_lower, xmax = IQR_higher, linetype = "IQR", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  # Add error bars for Range
  geom_errorbarh(
    aes(xmin = range_lower, xmax = range_higher, linetype = "Range", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +  
  # Add points for incubation periods, colored by subtype and shaped by estimate type
  geom_point(aes(color = Subtype, shape = value), size = 5, na.rm = TRUE) +
  # Customize colors for subtypes
  scale_color_manual(
    values = subtype_colors,
    name = "Subtype"  # Legend title for subtype colors
  ) +
  # Customize linetypes
  scale_linetype_manual(
    values = interval_linetypes,
    name = "Interval Type",  # Legend title for interval types
    labels = c("95% CI", "IQR", "Range")
  ) +
  # Customize shapes for estimate types
  scale_shape_manual(
    values = estimate_shapes,
    name = "Estimate Type",  # Legend title for estimate types
    breaks = c("Mean", "Median")  # Exclude NA from the legend
  ) +
  # Adjust X-axis scale
  scale_x_continuous(
    limits = c(0, 12),  # Set X-axis from 0 to 12
    breaks = seq(0, 12, by = 1)  # Tick marks in increments of 1
  ) +
  # Customize labels and titles
  labs(
    title = "",
    x = "Incubation Period (Days since exposure)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),         # General text size
    axis.text.y = element_text(size = 16, face = "bold"), # Y-axis text size
    axis.text.x = element_text(size = 16, face = "bold"), # Bold X-axis text
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(1.5, "cm")      # Adjust the overall size of legend keys
  )

ggsave(
  filename = file.path("plots", "inc_review.png"),
  plot = inc,
  width = 14,      
  height = 12,     
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)



###############################################################

#CFR

# Load the dataset (replace with your actual file path)
data <- read_excel(file.path("data", "H5N1pptdat.xlsx"), sheet = "CFR")


# Clean and standardize column names
colnames(data) <- iconv(colnames(data), to = "ASCII//TRANSLIT")  # Remove non-ASCII characters
colnames(data) <- gsub("[^[:print:]]", "", colnames(data))       # Replace non-printable characters
colnames(data) <- trimws(colnames(data))                        # Remove leading and trailing spaces

# Prepare the dataset for plotting
filtered_data <- data %>%
  mutate(
    Plot_Label = paste(trimws(Author), "-", trimws(Outbreak)),  # Clean and combine Author and Outbreak
    Fatality_risk = as.numeric(Fatality_risk),                 # Ensure numeric format for plotting
    IQR_lower = as.numeric(IQR_lower),
    IQR_higher = as.numeric(IQR_higher),
    `95%_lower` = as.numeric(`95%_lower`),
    `95%_higher` = as.numeric(`95%_higher`),
    range_lower = as.numeric(range_lower),
    range_higher = as.numeric(range_higher),
    Deaths_Cases = ifelse(is.na(Deaths) | is.na(Cases), NA, paste0(Deaths, "/", Cases))  # Combine Deaths/Cases
  )

# Correct and define custom order
custom_order <- rev(c(
  "CDC 2024 - 2024/25 US H5N1 outbreak (as of April 25)",
  "Serological studies - 2024/25 US H5N1 outbreak (as of April 25)",
  "Lai et al. 2016 - All casess 1997-2015",
  "Lai et al. 2016 - Clade 0",
  "Lai et al. 2016 - Clade 1",
  "Lai et al. 2016 - Clade 2.1",
  "Lai et al. 2016 - Clade 2.2",
  "Lai et al. 2016 - Clade 2.3",
  "Lai et al. 2016 - Clade 7",
  "Oner et al. 2012 - 0-5yrs (1997-2010)",
  "Oner et al. 2012 - 6-11yrs (1997-2010)",
  "Oner et al. 2012 - 12-17yrs (1997-2010)",
  "Oner et al. 2012 - All Pediatric Cases (1997-2010)",
  "Li et al. 2008 - Turkey  2006",
  "Li et al. 2008 - 1997 Hong Kong",
  "Fouchier et al. 2004 - Netherlands 2003",
  "Bosman et al. 2004 - Netherlands 2003",
  "McDonald et al. 2023 - Netherlands (2011/2012 - 2019/2020 seasons) (H3N2)",
  "McDonald et al. 2023 - Netherlands (2011/2012 - 2019/2020 seasons) (pdm09)",
  "McDonald et al. 2023 - Netherlands (2011/2012 - 2019/2020 seasons) (FLUBV)",
  "WHO 2017 - 1918 Pandemic",
  "WHO 2017 - 1957 Pandemic",
  "WHO 2017 - 1968 Pandemic",
  "WHO 2017 - 2009 Pandemic",
  "Wong et al. 2013 - 2009 Pandemic - laboratory-confirmed cases",
  "Wong et al. 2013 - 2009 Pandemic -symptomatic cases",
  "Wong et al. 2013 - 2009 Pandemic -infections",
  "Wong et al. 2013 - 2009 Pandemic -symptomatic cases (children)",
  "Wong et al. 2013 - 2009 Pandemic -symptomatic cases (elderly)",
  "Abdalla et al. 2020 - 2009 Pandemic (Saudi Arabia - hospitalised cases)",
  "Fujikura et al. 2014 - 2009 Pandemic (adult patients with pneumonia - Japan)",
  "Petrovic et al. 2011 - 2009 Pandemic (hospitalised cases - Serbia)", 
  "Presanis et al. 2009 - 2009 Pandemic (symptomatic cases - US)",
  "Presanis et al. 2009 - 2009 Pandemic (ICU cases - US)",
  "Presanis et al. 2009 - 2009 Pandemic (hospitalised cases - US)",
  "Larrauri Camara et al. 2012 - 2009 Pandemic (Spain)",
  "Altmann et al. 2011 - 2009 Pandemic (pediatric intensive care - Germany)",
  "Kim et al. 2023 - Influenza seasons  2010 to  2017 (hospitalised patients - Canada)"
))

# Apply the custom reversed order to the Plot_Label column
filtered_data <- filtered_data %>%
  mutate(Plot_Label = factor(Plot_Label, levels = custom_order))

# Define subtype colors
subtype_colors <- c(
  "H1N1" = "darkgreen", 
  "H2N2" = "darkblue", 
  "H3N2" = "orange", 
  "H1N1pdm09" = "grey", 
  "H5N1" = "red", 
  "Seasonal (mixed subtypes)" = "magenta",
  "H7N7" = "#f26b20",
  "Influenza B" = "brown"
)

# Create the forest plot
CFR <- ggplot(filtered_data, aes(x = Fatality_risk, y = Plot_Label)) +
  # Add error bars for IQR with subtype-colored lines
  geom_errorbarh(
    aes(xmin = IQR_lower, xmax = IQR_higher, linetype = "IQR", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  # Add error bars for 95% CI with subtype-colored lines
  geom_errorbarh(
    aes(xmin = `95%_lower`, xmax = `95%_higher`, linetype = "95% CI", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE, alpha = 0.6
  ) +
  # Add error bars for Range with subtype-colored lines
  geom_errorbarh(
    aes(xmin = range_lower, xmax = range_higher, linetype = "Range", color = Subtype),
    height = 0.3, size = 1.5, na.rm = TRUE, alpha = 0.8
  ) +
  # Add points for fatality risk
  geom_point(aes(color = Subtype, shape = Risk_type), size = 5, na.rm = TRUE) +
  # Add deaths/cases as text annotation
  #geom_text(
  # data = filtered_data %>% filter(!is.na(Deaths_Cases)),
  #aes(label = Deaths_Cases),
  #hjust = -0.5, vjust = -0.5, size = 3.5, na.rm = TRUE, color = "black"
  #) +
  # Customize colors for subtypes
  scale_color_manual(values = subtype_colors) +
  # Customize shapes for risk types
  scale_shape_manual(values = c("CFR" = 16, "IFR" = 17)) +  # Circle for CFR, Triangle for IFR
  # Customize linetypes for intervals
  scale_linetype_manual(
    values = c("IQR" = "solid", "95% CI" = "longdash", "Range" = "solid"),
    guide = guide_legend(title = "Interval Type")
  ) +
  # Customize labels and titles
  labs(
    title = "",
    x = "Fatality Risk (%)",
    y = "",
    color = "Subtype",
    shape = "Risk Type"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),         # General text size
    axis.text.y = element_text(size = 14), # Y-axis text size
    axis.text.x = element_text(size = 16), # Bold X-axis text
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(1.5, "cm")      # Adjust the overall size of legend keys
  )

ggsave(
  filename = file.path("plots", "CFR_review.png"),
  plot = CFR,
  width = 16,      
  height = 12,     
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)
#########################################################

#Serology 
data <- read_excel(file.path("data", "H5N1pptdat.xlsx"), sheet = "Sero")

# Clean column names and prepare the data
colnames(data) <- trimws(colnames(data))

data <- data %>%
  mutate(
    Plot_Label = paste(Author, "-", Population),  # Create labels for the y-axis
    Positive_Cases = ifelse(!is.na(Positive) & !is.na(Cases), paste0(Positive, "/", Cases), ""),  # Combine Positive and Cases
    Seroprevalence = as.numeric(Seroprevalence),  # Ensure numeric format
    `95%_lower` = as.numeric(`95%_lower`), 
    `95%_higher` = as.numeric(`95%_higher`)
  ) %>%
  arrange(Seroprevalence)  # Sort by Seroprevalence for better visualization

# Define custom colors and shapes
label_colors <- c(
  "1997-2020" = "black",
  "2024 US outbreak" = "red",
  "Egypt 2022-2023" = "blue"
)

data <- data %>%
  mutate(criteria = recode(criteria, "NS" = "Non-standard"))

criteria_shapes <- c("WHO" = 16, "Non-standard" = 17)  # Circle for WHO, Triangle for NS
custom_order <- c(
  "Chen et al. 2020 - General population (NS)",
  "Chen et al. 2020 - General population",
  "Chen et al. 2020 - Exposed health care workers (NS)",
  "Chen et al. 2020 - Exposed health care workers",
  "Chen et al. 2020 - Social contacts (NS)",
  "Chen et al. 2020 - Social contacts",
  "Chen et al. 2020 - Household contacts (NS)",
  "Chen et al. 2020 - Household contacts",
  "Chen et al. 2020 - Poultry cullers (NS)",
  "Chen et al. 2020 - Poultry cullers",
  "Chen et al. 2020 - Poultry workers (NS)",
  "Chen et al. 2020 - Poultry workers",
  "Gomaa et al. 2023 - Live bird market workers ( Egypt 2022-2023)",
  "Shittu et al. 2024 - Texas dairy farm workers (2024)",
  "Mellis et al. 2024 - Michigan and Colorado dairy workers (2024)"
)

data <- data %>%
  mutate(Plot_Label = factor(Plot_Label, levels = custom_order))



# Plot with the updated Y-axis order
sero <- ggplot(data, aes(x = Seroprevalence, y = Plot_Label)) +
  # Add error bars for 95% confidence intervals
  geom_errorbarh(
    aes(xmin = `95%_lower`, xmax = `95%_higher`, linetype = "95% CI"),
    height = 0.3, size = 1.2, color = "black", na.rm = TRUE
  ) +
  # Add points for seroprevalence with color based on Label and shape based on criteria
  geom_point(
    aes(color = Label, shape = criteria), 
    size = 5, 
    na.rm = TRUE 
  ) +
  # Customize colors and shapes
  scale_color_manual(
    values = label_colors,
    name = ""  # Legend title for color
  ) +
  scale_shape_manual(
    values = criteria_shapes,
    name = "Seropositivity criteria"  # Legend title for shape
  ) +
  scale_linetype_manual(
    values = c("95% CI" = "longdash"),
    name = "Interval Type"  # Title for line legend
  ) +
  # Customize labels and titles
  labs(
    title = "",
    x = "Seroprevalence (%)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),         # General text size
    axis.text.y = element_text(size = 16), # Y-axis text size
    axis.text.x = element_text(size = 16), # Bold X-axis text
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(1.5, "cm")      # Adjust the overall size of legend keys
  )

ggsave(
  filename = file.path("plots", "sero_review.png"),
  plot = sero,
  width = 14,      
  height = 12,     
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)


#########################################################
# infectious period

# Load the dataset (replace with your actual file path)
data <- read_excel(file.path("data", "H5N1pptdat.xlsx"), sheet = "Inf")

# Clean column names to remove unwanted spaces or characters
colnames(data) <- trimws(colnames(data))

# Prepare the data for plotting
data <- data %>%
  mutate(
    Plot_Label = paste(Author, "-", Outbreak),  # Create labels for the y-axis
    Mean_infectious_period = as.numeric(Mean_infectious_period),  # Ensure numeric format
    `95%_lower` = as.numeric(`95%_lower`), 
    `95%_higher` = as.numeric(`95%_higher`),
    range_lower = as.numeric(range_lower),
    range_higher = as.numeric(range_higher)
  ) %>%
  arrange(Mean_infectious_period)  # Sort by Mean infectious period for better visualization

# Define subtype colors
subtype_colors <- c(
  "H1N1/H3N2" = "brown", 
  "H3N2" = "orange", 
  "H5N1" = "red", 
  "Seasonal (mixed subtypes)" = "magenta",
  "H1N1pdm09" = "grey",
  "H1N1" = "darkgreen"
)

custom_order <- c(
  "Cauchemez et al. 2004 - 1999–2000 Influenza season (France)",
  "Tuite et al. 2010 - 2009 Pandemic (Canada)",
  "Lee et al. 2021 - 2009 Pandemic (Korea)",
  "Canini and Carratt 2011 - Volunteer challenge (H1N1)",
  "Cori et al. 2012 - Volunteer challenge",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Long inc.",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Primary inc.",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Short inc.",
  "Yang et al. 2007 - Indonesia 2006 (H5N1)"
)


# Apply the custom order to the data
data <- data %>%
  mutate(
    Plot_Label = factor(Plot_Label, levels = custom_order)  # Apply custom order
  )

inf <- ggplot(data, aes(x = Mean_infectious_period, y = Plot_Label)) +
  # Add error bars for 95% confidence intervals
  geom_errorbarh(
    aes(xmin = `95%_lower`, xmax = `95%_higher`, linetype = "95% CI", color = Subtype),
    height = 0.2, size = 1.2, na.rm = TRUE
  ) +
  # Add error bars for range
  geom_errorbarh(
    aes(xmin = range_lower, xmax = range_higher, linetype = "Range", color = Subtype),
    height = 0.2, size = 1.5, na.rm = TRUE, alpha = 0.6
  ) +
  # Add points for mean infectious period
  geom_point(aes(color = Subtype), size = 4, na.rm = TRUE) +
  # Customize colors for subtypes
  scale_color_manual(
    values = subtype_colors,
    name = "Subtype"  # Legend title for subtype colors
  ) +
  # Customize linetypes
  scale_linetype_manual(
    values = c("95% CI" = "longdash", "Range" = "solid"),
    name = "Line Type"  # Title for the line type legend
  ) +
  # Customize labels and titles
  labs(
    title = "Infectious Period",
    x = "Mean Infectious Period (Days)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),         # General text size
    axis.text.y = element_text(size = 13),                # Y-axis text size
    axis.text.x = element_text(size = 16),                # X-axis text size
    legend.text = element_text(size = 16, face = "bold"), # Legend text size
    legend.title = element_text(size = 18, face = "bold"),# Legend title size
    legend.key.size = unit(1.5, "cm")                     # Legend key size
  )

# Print plot
inf

###########################################################################################
#Latent period

data <- read_excel(file.path("data", "H5N1pptdat.xlsx"), sheet = "Lat")

# Clean column names to remove unwanted spaces or characters
colnames(data) <- trimws(colnames(data))

# Prepare the data for plotting
data <- data %>%
  mutate(
    Plot_Label = paste(Author, "-", Outbreak),  # Create labels for the y-axis
    Mean_Latent_period = as.numeric(Mean_Latent_period),  # Ensure numeric format
    `95%_lower` = as.numeric(`95%_lower`), 
    `95%_higher` = as.numeric(`95%_higher`),
    range_lower = as.numeric(range_lower),
    range_higher = as.numeric(range_higher)
  ) %>%
  arrange(Mean_Latent_period)  # Sort by Mean latent period for better visualization

# Define subtype colors
subtype_colors <- c(
  "H1N1/H3N2" = "brown", 
  "H1N1pdm09" = "grey", 
  "H5N1" = "red", 
  "Seasonal (mixed subtypes)" = "magenta",
  "H1N1" = "darkgreen"
)


custom_order <- c(
  "Tuite et al. 2010 - 2009 Pandemic (Canada)",
  "Canini and Carratt 2011 - Volunteer challenge (H1N1)",
  "Cori et al. 2012 - Volunteer challenge",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Long inc.",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Primary inc.",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Short inc."
)

# Apply the custom order to the data
data <- data %>%
  mutate(
    Plot_Label = factor(Plot_Label, levels = custom_order)  # Apply custom order
  )

# Create the latent period plot
lat <- ggplot(data, aes(x = Mean_Latent_period, y = Plot_Label)) +
  # Add error bars for 95% confidence intervals, colored by subtype
  geom_errorbarh(
    aes(xmin = `95%_lower`, xmax = `95%_higher`, linetype = "95% CI", color = Subtype),
    height = 0.2, size = 1.2, na.rm = TRUE
  ) +
  # Add error bars for range, colored by subtype
  geom_errorbarh(
    aes(xmin = range_lower, xmax = range_higher, linetype = "Range", color = Subtype),
    height = 0.2, size = 1.5, na.rm = TRUE, alpha = 0.6
  ) +
  # Add points for mean latent period colored by subtype
  geom_point(aes(color = Subtype), size = 4, na.rm = TRUE, show.legend = FALSE) + # Suppress color legend
  # Customize colors for subtypes
  scale_color_manual(
    values = subtype_colors,
    name = "Subtype"  # Legend title for subtype colors
  ) +
  # Customize linetypes
  scale_linetype_manual(
    values = c("95% CI" = "longdash", "Range" = "solid"),
    name = "Line Type",  # Title for the line type legend
    labels = c("95% CI", "Range")
  ) +
  # Customize labels and titles
  labs(
    title = "Latent Period",
    x = "Mean Latent Period (Days)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),         # General text size
    axis.text.y = element_text(size = 13),                # Y-axis text size
    axis.text.x = element_text(size = 16),                # X-axis text size
    legend.text = element_text(size = 16, face = "bold"), # Legend text size
    legend.title = element_text(size = 18, face = "bold"),# Legend title size
    legend.key.size = unit(1.5, "cm")                     # Legend key size
  ) +
  guides(color = "none")  # Suppress the color legend
lat

#combine plot
combined_plot <- lat / inf + plot_layout(guides = "collect")

# Print the combined plot
print(combined_plot)

ggsave(
  filename = file.path("plots", "lat_inf_review.png"),
  plot = combined_plot,
  width = 14,      
  height = 12,     
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

######################################################################

#Serial interval 
data <- read_excel(file.path("data", "H5N1pptdat.xlsx"), sheet = "SI")


subtype_colors <- c(
  "H1N1pdm09" = "grey",
  "H1N1" = "darkgreen", 
  "H3N2" = "orange",
  "H5N1" = "red",
  "Influenza B" = "blue",
  "Seasonal (mixed subtypes)" = "magenta",
  "H7N9" = "#ff69b4"
)


# Clean column names to remove unwanted spaces or characters
colnames(data) <- trimws(colnames(data))

# Define custom order for the Y-axis
custom_order <- c(
  "Vink et al. 2014 - Influenza B",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Long inc",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Primary  inc",
  "Chan et al. 2024 - 2021-2023 Influenza seasons (USA) - Short inc",
  "Vink et al. 2014 - 2009 Pandemic (Netherlands)",
  "Vink et al. 2014 - 2009 Pandemic (Canada) (2)",
  "Vink et al. 2014 - 2009 Pandemic (Canada) (1)",
  "Vink et al. 2014 - 2009 Pandemic (US) (3)",
  "Vink et al. 2014 - 2009 Pandemic (US) (2)",
  "Vink et al. 2014 - 2009 Pandemic (US) (1)",
  "Vink et al. 2014 - H1N1pdm09",
  "Vink et al. 2014 - 1999–2000 Influenza season (France)",
  "Vink et al. 2014 - H3N2",
  "Canini and Carratt 2011 - Volunteer challenge",
  "Vink et al. 2014 - H1N1",
  "Estimated - Indonesia 2005–2009 (Log normal)",
  "Estimated - Indonesia 2005–2009 (Gamma)"
)

# Prepare the data for plotting
data <- data %>%
  mutate(
    Plot_Label = paste(Author, "-", Outbreak),  # Remove subtype and estimate type from label
    Plot_Label = factor(Plot_Label, levels = custom_order),  # Apply custom order
    Serial_interval = as.numeric(Serial_interval),  # Ensure numeric format
    `95%_CI_lower` = as.numeric(`95%_CI_lower`),
    `95%_CI_higher` = as.numeric(`95%_CI_higher`),
    `95%_CrI_lower` = as.numeric(`95%_CrI_lower`),
    `95%_CrI_higher` = as.numeric(`95%_CrI_higher`),
    `range_lower` = as.numeric(`range_lower`),
    `range_higher` = as.numeric(`range_higher`)
  )



# Define shapes for estimate types
estimate_shapes <- c("Mean" = 17, "Median" = 16)  # Triangle for Mean, Circle for Median

# Define interval linetypes
interval_linetypes <- c(
  "95% CI" = "longdash",
  "95% CrI" = "dashed",
  "Range" = "solid"
)

# Dynamically adjust the X-axis limits to include all data points
x_limits <- range(c(
  data$`95%_CI_lower`, data$`95%_CI_higher`,
  data$`95%_CrI_lower`, data$`95%_CrI_higher`,
  data$`range_lower`, data$`range_higher`
), na.rm = TRUE)

serial <- ggplot(data, aes(x = Serial_interval, y = Plot_Label)) +
  # Error bars (same as before)
  geom_errorbarh(
    aes(xmin = `95%_CI_lower`, xmax = `95%_CI_higher`, linetype = "95% CI", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  geom_errorbarh(
    aes(xmin = `range_lower`, xmax = `range_higher`, linetype = "Range", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  geom_errorbarh(
    aes(xmin = `95%_CrI_lower`, xmax = `95%_CrI_higher`, linetype = "95% CrI", color = Subtype),
    height = 0.3, size = 1.2, na.rm = TRUE
  ) +
  # Plot points, but drop NAs from shape legend with show.legend = FALSE for those
  geom_point(
    data = data[!is.na(data$Estimate_type), ],
    aes(color = Subtype, shape = Estimate_type),
    size = 5
  ) +
  geom_point(
    data = data[is.na(data$Estimate_type), ],
    aes(color = Subtype),
    shape = 16, size = 5, show.legend = FALSE  # use circle or any shape
  ) +
  # Customize colors, shapes, linetypes (same as before)
  scale_color_manual(values = subtype_colors, name = "Subtype") +
  scale_shape_manual(values = estimate_shapes, name = "Estimate Type", na.translate = FALSE) +
  scale_linetype_manual(
    values = interval_linetypes,
    name = "Interval Type",
    labels = c("95% CI", "95% CrI", "Range"),
    guide = guide_legend(order = 1)
  ) +
  guides(
    shape = guide_legend(order = 2),
    color = guide_legend(order = 3)
  ) +
  coord_cartesian(xlim = c(0, 10)) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(
    title = "",
    x = "Serial Interval (Days)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 16),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(1.5, "cm")
  )


ggsave(
  filename = file.path("plots", "serial_review.png"),
  plot = serial,
  width = 14,      
  height = 12,     
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)


#######################################################################################