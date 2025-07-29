library(shiny)

# dev.off() # Sometimes need this to get the shiny to show the pheatmap?


source("Import_data.R") # for list_dfs_2 and my_tpm
source("Import_GeneSets.R")

# Adjust the dataframes
my_tpm <- All_tpm
names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it
my_tpm <- my_tpm %>% column_to_rownames("X")

pipeSummary_2 <- All_pipeSummary %>% filter(SampleID %in% c(SputumSampleList, BrothSampleList, MarmSampleList, MimicSampleList, RabbitSampleList)) %>%mutate(SampleID = gsub("_S.*", "", SampleID))


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=12), 
        plot.margin = margin(10, 10, 10, 20))

my_annotation_colors <- list(
  Type2 = c("W0 sputum (cure)" = "#0072B2",
            "W0 sputum (relapse)" = "red", 
            "Caseum mimic" = "green4",
            "Marmoset" = "#6A3D9A", 
            "Rabbit" = "#E69F00", 
            "Broth" = "#999999")
)


# Define UI ----
ui <- fluidPage(
  titlePanel("Precious Sputum Run 1 Pheatmap"),
  
  fluidRow(
    
    column(width = 4,
           
           # Which sample types to plot
           checkboxGroupInput("SampleTypes",
                              label = "Select sample type",
                              choices = c("W0 sputum (cure)", "W0 sputum (relapse)", "Caseum mimic", "Marmoset", "Rabbit", "Broth"),
                              selected = c("W0 sputum (cure)", "W0 sputum (relapse)", "Caseum mimic", "Marmoset", "Rabbit", "Broth")),
           
           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList)),
           # Dropdown for selecting the gene set within the chosen rda file.
           selectInput("my_GeneSet",
                       label = "Gene Set",
                       choices = NULL),
           
           # Select genes manually
           textInput("manual_genes",
                     label = "Enter Rv# (comma separated)",
                     placeholder = "Rv1473A, Rv2011c, Rv0494"),
           # Add row clustering options
           numericInput("cutree_rows", 
                        label = "Number of Row Clusters", 
                        value = 1, min = 1, step = 1),
           # Add column clustering options
           numericInput("cutree_cols", 
                        label = "Number of Column Clusters", 
                        value = 1, min = 1, step = 1),
           selectInput("my_scaling",
                       label = "How to scale",
                       choices = (c("row", "column", "none"))),
           # Add checkbox to toggle display_numbers in heatmap
           checkboxInput("show_numbers", label = "Show Values", value = FALSE),
           textInput("my_GeneID", 
                     label = "Pick a gene to link to mycobrowser",
                     placeholder = "Rv..."),
           uiOutput("gene_link")  # New UI output for the link
    ),
    
    column(width = 8, # Max is 12...
           uiOutput("dynamic_pheatmap")
           # plotOutput("pheatmap", width = "200%", height = "600")
           # plotOutput("pheatmap", width = "100%", height = "600px")
    )
    
  )
  
)

# Define server logic ----
server <- function(input, output, session) {
  
  # Gene Link
  output$gene_link <- renderUI({
    req(input$my_GeneID)  # Ensure there's a valid input
    url <- paste0("https://mycobrowser.epfl.ch/genes/", input$my_GeneID)
    tags$a(href = url, target = "_blank", paste0("View Details of ", input$my_GeneID, " on Mycobrowser"))
  })
  
  # When a new gene set source is selected, update the gene set dropdown
  observeEvent(input$my_GeneSetSource, {
    updateSelectInput(session, "my_GeneSet",
                      choices = names(allGeneSetList[[input$my_GeneSetSource]]),
                      selected = NULL)
  })
  
  # Function to get the selected genes
  get_selected_genes <- reactive({
    if (input$manual_genes != "") {
      # Process manual input: remove extra spaces, split by commas or spaces
      genes <- unlist(strsplit(input$manual_genes, "[,\\s]+"))
      genes <- trimws(genes)  # Trim whitespace
      genes <- genes[genes != ""]  # Remove empty entries
    } else if (!is.null(input$my_GeneSet) && input$my_GeneSet %in% names(allGeneSetList[[input$my_GeneSetSource]])) {
      genes <- allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]]
    } else {
      genes <- character(0)  # Empty vector if nothing is selected
    }
    
    return(genes)
  })
  
  # Dynamic UI for heatmap height
  output$dynamic_pheatmap <- renderUI({
    # req(input$my_GeneSetSource, input$my_GeneSet)
    
    # Get the list of genes from either dropdown selection or manual input
    selected_genes <- get_selected_genes()
    
    # Count the number of genes in the selected set
    num_genes <- sum(rownames(my_tpm) %in% allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]])
    
    # Dynamically set plot height (base height + extra space per gene)
    plot_height <- max(400, min(2500, num_genes * 40))  # Adjust as needed
    # plot_width <- max(2000, min(2000)) # Still need to figure out how to change this! 
    
    plotOutput("pheatmap", height = paste0(plot_height, "px"))
  })
  
  
  # Render the pheatmap
  output$pheatmap <- renderPlot({
    
    # Make a list of the columns for each timepoint
    type_columns <- list(
      "W0 sputum (cure)" = c("W0_11011", "W0_11012", "W0_11047", "W0_12007", "W0_12008", "W0_12009", "W0_12010", "W0_12012", "W0_12019", "W0_12028", "W0_12031", "W0_12032", "W0_12034", "W0_12070", "W0_12083", "W0_13004", "W0_13012", "W0_13016", "W0_13045", "W0_14005", "W0_14063", "W0_14113", "W0_15017", "W0_15018", "W0_15029", "W0_15045", "W0_15065", "W0_15081", "W0_15089"),
      "W0 sputum (relapse)" = c("W0_11031", "W0_11058", "W0_12073", "W0_13041", "W0_14025"),
      "Caseum mimic" = c("HN878_mimic_D14_R1", "HN878_mimic_D14_R2", "HN878_mimic_D28_R1", "HN878_mimic_D7_R2"),
      "Marmoset" = c("BQ12_10_Probe_3A", "BQ12_8_Probe_4A_50"),
      "Rabbit" = c("LLL_Cav_L2_18", "LLL_Cav_L2_21", "LU_Cav_L1_12"),
      "Broth" = c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6"))
    
    selected_genes <- get_selected_genes()
    # print(selected_genes)
    
    # Filter data based on selected genes
    my_data <- my_tpm[rownames(my_tpm) %in% selected_genes, , drop = FALSE]
    
    # Now filter columns based on which SampleTypes are checked
    if (!is.null(input$SampleTypes)) {
      selected_cols <- unlist(type_columns[input$SampleTypes])
      selected_cols <- intersect(selected_cols, colnames(my_data))  # only keep existing columns
      my_data <- my_data[, selected_cols, drop = FALSE]
    }
    
    # If no columns left, show error
    if (ncol(my_data) == 0) {
      showNotification("No columns match the selected SampleTypes.", type = "error")
      return(NULL)
    }
    
    # Check if we have at least 2 genes
    if (nrow(my_data) < 2) {
      showNotification("At least two valid genes are required for clustering.", type = "error")
      return(NULL)
    }
    
    
    # Sort out the color categories
    my_color_annotation <- pipeSummary_2 %>%
      filter(SampleID %in% colnames(my_data)) %>%
      column_to_rownames("SampleID") %>%
      select(Type2)
    
    
    p <- pheatmap(my_data, 
                  annotation_col = my_color_annotation, 
                  annotation_colors = my_annotation_colors,
                  # col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                  scale = input$my_scaling, 
                  display_numbers = input$show_numbers,
                  fontsize_number = 8,
                  cutree_rows = input$cutree_rows,
                  cutree_cols = input$cutree_cols,
                  fontsize = 18)
    p
    
    # pheatmap returns a complex grid object; use grid.draw() to render it in Shiny.
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)# , options = list(launch.browser = TRUE))