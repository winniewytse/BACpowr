#' An Interactive Shiny App for HCB Power Analysis
#'
#'
#'
hcbApp <- function() {
  #### UI ####

  ui <- fluidPage(
    navbarPage("Hybrid Classical-Bayesian Power Analysis",
               navbarMenu("Study Design",
                          tabPanel("2-Level Cluster Randomized",
                                   sidebarLayout(
                                     sidebarPanel(
                                       tags$style(".well {background-color: #EEFEF2;}"),
                                       h3("Required Inputs"),
                                       numericInput("d_est", h5("Effect Size Estimate"),
                                                    value = .35),
                                       numericInput("d_sd", h5("SE of Effect Size Estimate"),
                                                    value = 0, min = 0),
                                       numericInput("rho_est", h5("ICC Estimate"),
                                                    value = .15, min = 0, max = 1),
                                       numericInput("rho_sd", h5("SE of ICC Estimate"),
                                                    value = 0, min = 0, max = 1),
                                       numericInput("r2_est", h5("Cluster-Level Covariate R2 Estimate"),
                                                    value = 0, min = 0, max = 1),
                                       numericInput("r2_sd", h5("SE of Cluster-Level Covariate R2 Estimate"),
                                                    value = 0, min = 0, max = 1),
                                       # numericInput("n", h5("Cluster Size"),
                                       #              value = NULL, min = 1),
                                       numericInput("K", h5("Number of Cluster-Level Covariate"),
                                                    value = 0, min = 0),
                                       sliderInput("power", h5("Desired Power Level"),
                                                   value = .80, min = 0, max = 1),
                                       h3("Optional Inputs"),
                                       selectInput("fixedJn", "Is J or n fixed?",
                                                   choices = c("No",
                                                               "Fixed number of clusters (J)",
                                                               "Fixed cluster size (n)"),
                                                   selected = "Fixed cluster size (n)"),
                                       conditionalPanel(
                                         condition = "input.fixedJn == 'Fixed number of clusters (J)'",
                                         numericInput("J", h5("J"),
                                                      value = 50, min = 3)
                                       ),
                                       conditionalPanel(
                                         condition = "input.fixedJn == 'Fixed cluster size (n)'",
                                         numericInput("n", h5("n"),
                                                      value = 20, min = 3)
                                       )
                                     ),
                                     mainPanel(
                                       plotOutput("plot1", height = "800px") %>%
                                         shinycssloaders::withSpinner(type = 6, size = .8, color = "#1DD6A2"),
                                       br(),
                                       br(),
                                       textOutput("est"),
                                       br(),
                                       textOutput("Jn")
                                     )
                                   )
                          )
               )
    )
  )

  #### SERVER ####

  server <- function(input, output) {

    output$plot1 <- output$plot2 <- output$plot3 <- renderPlot({
      if (input$fixedJn == "Fixed number of clusters (J)") {
        crtJn(d_est = input$d_est, d_sd = input$d_sd,
              rho_est = input$rho_est, rho_sd = input$rho_sd,
              r2_est = input$r2_est, r2_sd = input$r2_sd,
              J = input$J,K = input$K,
              power = input$power, plot = TRUE)[[1]]
      } else if (input$fixedJn == "Fixed cluster size (n)") {
        crtJn(d_est = input$d_est, d_sd = input$d_sd,
              rho_est = input$rho_est, rho_sd = input$rho_sd,
              r2_est = input$r2_est, r2_sd = input$r2_sd,
              n = input$n,K = input$K,
              power = input$power, plot = TRUE)[[1]]
      } else {
        crtJn(d_est = input$d_est, d_sd = input$d_sd,
              rho_est = input$rho_est, rho_sd = input$rho_sd,
              r2_est = input$r2_est, r2_sd = input$r2_sd,
              K = input$K,
              power = input$power, plot = TRUE)[[1]]
      }
    })

    output$est <- renderText({
      paste0(
        "You provided an effect size estimate of ", input$d_est,
        " with a SE of ", input$d_sd,
        ", an intraclass correlation estaimte of ", input$rho_est,
        " with a SE of ", input$rho_sd,
        ", and an estimate of variances explained by ", input$K,
        " cluster-level covariate(s) of ",
        input$r2_est,
        " with a SE of ", input$r2_sd, ". "
      )
    })

    output$Jn <- renderText({
      if (input$fixedJn == "Fixed number of clusters (J)") {
        Jn <- crtJn(d_est = input$d_est, d_sd = input$d_sd,
                    rho_est = input$rho_est, rho_sd = input$rho_sd,
                    r2_est = input$r2_est, r2_sd = input$r2_sd,
                    J = input$J, K = input$K,
                    power = input$power, plot = TRUE)
        paste0("For the given inputs, it requires J = ", Jn[[2]][1],
               " and n = ", Jn[[2]][2], " to achieve ",
               input$power*100, "% expected power. ")
      } else if (input$fixedJn == "Fixed cluster size (n)") {
        Jn <- crtJn(d_est = input$d_est, d_sd = input$d_sd,
                    rho_est = input$rho_est, rho_sd = input$rho_sd,
                    r2_est = input$r2_est, r2_sd = input$r2_sd,
                    n = input$n,K = input$K,
                    power = input$power, plot = TRUE)
        paste0("For the given inputs, it requires J = ", Jn[[2]][1],
               " and n = ", Jn[[2]][2], " to achieve ",
               input$power*100, "% expected power. ")
      } else {
        Jn <- crtJn(d_est = input$d_est, d_sd = input$d_sd,
                    rho_est = input$rho_est, rho_sd = input$rho_sd,
                    r2_est = input$r2_est, r2_sd = input$r2_sd,
                    # J = input$J, n = input$n,
                    K = input$K,
                    power = input$power, plot = TRUE)
        paste0("For the given inputs, it requires J = ", Jn[[2]][1],
               " and n = ", Jn[[2]][2], " to achieve ",
               input$power*100, "% expected power. ")
      }
    })
  }

  shinyApp(ui, server)
}

