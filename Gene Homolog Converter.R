# Load necessary libraries
library(shiny)
library(biomaRt)
library(DT) # For interactive tables

# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("\n      body {\n        background-color: lavender;\n        color: #4B0082;\n        font-family: Arial, sans-serif;\n      }\n      .shiny-input-container {\n        color: #4B0082;\n      }\n      .dataTables_wrapper .dataTables_paginate .paginate_button {\n        color: #4B0082 !important;\n      }\n      .dataTables_wrapper .dataTables_length select,\n      .dataTables_wrapper .dataTables_filter input {\n        background-color: lavender;\n        color: #4B0082;\n        border: 1px solid #4B0082;\n      }\n      .btn-primary {\n        background-color: #4B0082;\n        border-color: #4B0082;\n        color: white;\n      }\n      .btn-primary:hover {\n        background-color: #6A5ACD;\n        border-color: #6A5ACD;\n      }\n    "))
  ),
  titlePanel("Gene Homolog Finder"),
  sidebarLayout(
    sidebarPanel(
      textInput("human_gene", "Enter Human Gene(s) (comma-separated)", value = ""),
      textInput("mouse_gene", "Enter Mouse Gene(s) (comma-separated)", value = ""),
      textInput("rat_gene", "Enter Rat Gene(s) (comma-separated)", value = ""),
      actionButton("find_homologs", "Find Homologs", class = "btn-primary")
    ),
    mainPanel(
      h3("Results"),
      DTOutput("homolog_table")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  # BioMart connections
  mart_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mart_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mart_rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  
  # Reactive values for results
  homologs <- reactiveVal(data.frame())
  
  # Action on button click
  observeEvent(input$find_homologs, {
    message("Button clicked!")
    
    # Split comma-separated inputs
    human_genes <- strsplit(input$human_gene, ",")[[1]]
    mouse_genes <- strsplit(input$mouse_gene, ",")[[1]]
    rat_genes <- strsplit(input$rat_gene, ",")[[1]]
    
    # Trim whitespace
    human_genes <- trimws(human_genes)
    mouse_genes <- trimws(mouse_genes)
    rat_genes <- trimws(rat_genes)
    
    # Initialize results
    results <- data.frame(Human = NA, Mouse = NA, Rat = NA, stringsAsFactors = FALSE)
    
    tryCatch({
      if (length(human_genes) > 0 && nzchar(human_genes[1])) {
        for (human_gene in human_genes) {
          message("Querying for human gene: ", human_gene)
          human_gene_info <- getBM(
            attributes = c("hgnc_symbol", "ensembl_gene_id"),
            filters = "hgnc_symbol",
            values = human_gene,
            mart = mart_human
          )
          mouse_homologs <- getBM(
            attributes = c("ensembl_gene_id", "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_confidence"),
            filters = "ensembl_gene_id",
            values = human_gene_info$ensembl_gene_id,
            mart = mart_human
          )
          rat_homologs <- getBM(
            attributes = c("ensembl_gene_id", "rnorvegicus_homolog_associated_gene_name", "rnorvegicus_homolog_orthology_confidence"),
            filters = "ensembl_gene_id",
            values = human_gene_info$ensembl_gene_id,
            mart = mart_human
          )
          results <- rbind(results, data.frame(
            Human = paste0("<span style='color:red;'>", human_gene, "</span>"),
            Mouse = ifelse(nrow(mouse_homologs) > 0, paste0(mouse_homologs$mmusculus_homolog_associated_gene_name, " (", mouse_homologs$mmusculus_homolog_orthology_confidence, ")"), "NA"),
            Rat = ifelse(nrow(rat_homologs) > 0, paste0(rat_homologs$rnorvegicus_homolog_associated_gene_name, " (", rat_homologs$rnorvegicus_homolog_orthology_confidence, ")"), "NA")
          ))
        }
      }
      
      if (length(mouse_genes) > 0 && nzchar(mouse_genes[1])) {
        for (mouse_gene in mouse_genes) {
          message("Querying for mouse gene: ", mouse_gene)
          mouse_gene_info <- getBM(
            attributes = c("mgi_symbol", "ensembl_gene_id"),
            filters = "mgi_symbol",
            values = mouse_gene,
            mart = mart_mouse
          )
          human_homologs <- getBM(
            attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
            filters = "ensembl_gene_id",
            values = mouse_gene_info$ensembl_gene_id,
            mart = mart_mouse
          )
          rat_homologs <- getBM(
            attributes = c("ensembl_gene_id", "rnorvegicus_homolog_associated_gene_name", "rnorvegicus_homolog_orthology_confidence"),
            filters = "ensembl_gene_id",
            values = mouse_gene_info$ensembl_gene_id,
            mart = mart_mouse
          )
          results <- rbind(results, data.frame(
            Human = ifelse(nrow(human_homologs) > 0, paste0(human_homologs$hsapiens_homolog_associated_gene_name, " (1)"), "NA"),
            Mouse = paste0("<span style='color:red;'>", mouse_gene, "</span>"),
            Rat = ifelse(nrow(rat_homologs) > 0, paste0(rat_homologs$rnorvegicus_homolog_associated_gene_name, " (", rat_homologs$rnorvegicus_homolog_orthology_confidence, ")"), "NA")
          ))
        }
      }
      
      if (length(rat_genes) > 0 && nzchar(rat_genes[1])) {
        for (rat_gene in rat_genes) {
          message("Querying for rat gene: ", rat_gene)
          rat_gene_info <- getBM(
            attributes = c("rgd_symbol", "ensembl_gene_id"),
            filters = "rgd_symbol",
            values = rat_gene,
            mart = mart_rat
          )
          human_homologs <- getBM(
            attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
            filters = "ensembl_gene_id",
            values = rat_gene_info$ensembl_gene_id,
            mart = mart_rat
          )
          mouse_homologs <- getBM(
            attributes = c("ensembl_gene_id", "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_confidence"),
            filters = "ensembl_gene_id",
            values = rat_gene_info$ensembl_gene_id,
            mart = mart_rat
          )
          results <- rbind(results, data.frame(
            Human = ifelse(nrow(human_homologs) > 0, paste0(human_homologs$hsapiens_homolog_associated_gene_name, " (1)"), "NA"),
            Mouse = ifelse(nrow(mouse_homologs) > 0, paste0(mouse_homologs$mmusculus_homolog_associated_gene_name, " (", mouse_homologs$mmusculus_homolog_orthology_confidence, ")"), "NA"),
            Rat = paste0("<span style='color:red;'>", rat_gene, "</span>")
          ))
        }
      }
    }, error = function(e) {
      message("Error during BioMart query: ", e$message)
    })
    
    # Assign results to reactive value
    homologs(results)
  })
  
  # Render results as a styled table
  output$homolog_table <- renderDT({
    result_table <- homologs()
    datatable(result_table, escape = FALSE, options = list(dom = 't'))
  })
}

# Run the application
shinyApp(ui = ui, server = server)

