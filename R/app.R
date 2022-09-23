# Load packages
library(shiny)
library(shinythemes)
library(dplyr)
library(dqshiny)

# Load data
utils::data("disgenet")
utils::data("interactome_subset")
disgenet <- disgenet %>%
  filter(score >= 0.6)

# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
  titlePanel("Drug Synergy"),
  sidebarLayout(
    sidebarPanel(
      # Select disease 1
      autocomplete_input(id = "diseaseName1", label = strong("Disease 1"),
                  options = unique(disgenet$diseaseName),
                  value = "Pulmonary Hypertension"),
      # Select disease 2
      autocomplete_input(id = "diseaseName2", label = strong("Disease 2"),
                  options = unique(disgenet$diseaseName),
                  value = "Pulmonary Hypertension"),
      # Select drug 1
      autocomplete_input(id = "drug1", label = strong("Drug 1"),
                  options = unique(interactome_subset$Protein_A),
                  value = "Hydrochlorothiazide"),
      # Select drug 2
      autocomplete_input(id = "drug2", label = strong("Drug 2"),
                  options = unique(interactome_subset$Protein_B),
                  value = "Diazoxide")
    ),

    # Output: Description, lineplot, and reference
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Network Graph",
          plotOutput(outputId = "netgraph", width = "1000px", height = "1000px")
        )
      )
    )
  )
)

# Define server function
server <- function(input, output) {

  # Subset data
  selected_type <- reactive({
    req(input$diseaseName1)
    req(input$diseaseName2)
    req(input$drug1)
    req(input$drug2)
    # validate(need(input$diseaseName1 != input$diseaseName2, "Error: Diseases must be different."))
    load_subnetwork(input$diseaseName1, input$diseaseName2, input$drug1, input$drug2)
  })

  output$netgraph <- renderPlot({
    plot(selected_type()$subnetwork,
         vertex.color = RColorBrewer::brewer.pal(length(unique(V(selected_type()$subnetwork)$source)), "Dark2")[as.numeric(as.factor(vertex_attr(selected_type()$subnetwork, "source")))])
    })

}

# Create Shiny object
shinyApp(ui = ui, server = server)
