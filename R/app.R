# Load packages
library(shiny)
library(shinythemes)
library(dplyr)
library(dqshiny)

# Load data
utils::data("disgenet")
utils::data("interactome_subset")

# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
  titlePanel("Drug Synergy"),
  sidebarLayout(
    sidebarPanel(

      # Select type of trend to plot
      selectInput(inputId = "diseaseType", label = strong("Disease type"),
                  choices = unique(disgenet$diseaseType),
                  selected = "disease"),
      # Select disease 1
      selectInput(inputId = "diseaseName1", label = strong("Disease 1"),
                  choices = unique(disgenet$diseaseName),
                  selected = "Pulmonary Hypertension"),
      # Select disease 2
      selectInput(inputId = "diseaseName2", label = strong("Disease 2"),
                  choices = unique(disgenet$diseaseName),
                  selected = "Liver carcinoma"),
      # Select drug 1
      autocomplete_input(id = "drug1", label = strong("Drug 1"),
                  options = unique(interactome_subset$Protein_A),
                  value = "Imatinib"),
      # Select drug 2
      autocomplete_input(id = "drug2", label = strong("Drug 2"),
                  options = unique(interactome_subset$Protein_B),
                  value = "Tandutinib")
    ),

    # Output: Description, lineplot, and reference
    mainPanel(
      tableOutput(outputId = "table")
    )
  )
)

# Define server function
server <- function(input, output) {

  # Subset data
  selected_type <- reactive({
    req(input$diseaseName1)
    req(input$diseaseName2)
    validate(need(input$diseaseName1 != input$diseaseName2, "Error: Diseases must be different."))
    disgenet %>%
      filter(
        diseaseType == input$diseaseType
        )
  })

  output$table <- renderTable({
    head(selected_type())
    })

}

# Create Shiny object
shinyApp(ui = ui, server = server)
