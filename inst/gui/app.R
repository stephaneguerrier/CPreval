#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(CPreval)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Conditional Prevalence"),
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(

            numericInput("R",
                         h4("R"),
                         value = 5),

            numericInput("R0",
                         h4("R0"),
                         value = 3),

            numericInput("n",
                         h4("n"),
                         value = 1544),

            numericInput("pi0",
                         h4("pi0"),
                         value = 1/758),

            numericInput("gamma",
                         h4("gamma"),
                         value = 0.05),

            numericInput("alpha",
                         h4("alpha"),
                         value = 0),

            numericInput("alpha0",
                         h4("alpha0"),
                         value = 0),

            numericInput("beta",
                         h4("beta"),
                         value = 0),

            numericInput("beta0",
                         h4("beta0"),
                         value = 0)

        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("MLE", verbatimTextOutput("mle")),
                tabPanel("Survey sample", verbatimTextOutput("survey")),
                tabPanel("Moment Estimator", verbatimTextOutput("moment"))

            )
        )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
    output$mle <- renderPrint({

        mle(R = input$R, R0 = input$R0,
            pi0 = input$pi0, n = input$n,
            alpha = input$alpha, beta = input$beta,
            alpha0 = input$alpha0, beta0 = input$beta0)

    })

    output$moment <- renderPrint({

        moment_estimator(R = input$R, R0 = input$R0,
                         pi0 = input$pi0, n = input$n,
                         alpha = input$alpha, beta = input$beta,
                         alpha0 = input$alpha0, beta0 = input$beta0)

    })

    output$survey <- renderPrint({

        survey_sample(R = input$R,
                         pi0 = input$pi0, n = input$n,
                         alpha = input$alpha, beta = input$beta)

    })

}

# Run the application
shinyApp(ui = ui, server = server)
