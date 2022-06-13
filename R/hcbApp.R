#' An Interactive Shiny App for HCB Power Analysis
#'
#' \code{hcbApp()} executes the interactive Shiny App for HCB power analysis.
#'
#' @return An interactive Shiny App to perform HCB power analysis.
#' @import shiny
#' @import shinydashboard
#' @export
#'
hcbApp <- function() {

  #### UI ####

  ui <- dashboardPage(
    dashboardHeader(title = "Sample Size Planning with the HCB Approach",
                    titleWidth = 500),
    dashboardSidebar(width = 150,
                     sidebarMenu(
                       menuItem(
                         "Instruction",
                         tabName = "instruction",
                         icon = icon("home")
                       ),
                       menuItem(
                         "Two-Level CRT",
                         tabName = "crt2",
                         icon = icon("users")
                       )
                     )),
    dashboardBody(
      tabItems(
        tabItem(
          tabName = "instruction",
          fluidPage(
            box(
              title = "Welcome!",
              status = "primary",
              width = 12,
              solidHeader = FALSE,
              p(
                "This application plans sample size for a new study ",
                "using the hybrid classical-Bayesian (HCB) approach. The HCB approach ",
                "incorporates uncertainty in parameter values to more realistically ",
                "determine the sample size requisite. More features and study designs ",
                "will be added in the near future. Please stay tuned :)"
              ),
              p(
                "If you have any questions or suggestions to this application or ",
                "approach, please email wingyeet@usc.edu. If you find this application ",
                "helpful, please cite the following article. "
              ),
              p(
                "Thank you for visiting, and we wish you all the best in your research!"
              )
            ),
            box(
              title = "Reference",
              status = "primary",
              width = 12,
              solidHeader = FALSE,
              p(
                "Tse, W. W. & Lai, M. H. C. (in press). ",
                "Incorporating Uncertainty in Power Analysis: A Hybrid",
                "Classical-Bayesian Approach for Designing Two-Level Cluster ",
                "Randomized Trials.",  em("Psychological Methods. ")
              )
            )
          )
        ),
        tabItem(
          tabName = "crt2",
          titlePanel("Design a Two-Level Cluster Randomized Trial"),
          fluidPage(
            style = 'padding:0px;',
            box(
              width = 12,
              style = 'margin:0px;',
              status = "primary",
              h3("Inputs"),
              fluidRow(
                column(
                  width = 4,
                  withMathJax(),
                  numericInput(
                    inputId = "d_est",
                    label = h5("\\(\\delta\\) Effect size value"),
                    value = .5
                  ),
                  numericInput(
                    inputId = "rho_est",
                    label = h5("\\(\\rho\\) ICC value"),
                    value = .2,
                    min = 0, max = 1
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "d_sd",
                    label = h5("\\(\\sigma_\\delta\\)
                               Uncertatinty level of \\(\\delta\\)"),
                    value = .1, min = 0
                  ),
                  numericInput(
                    inputId = "rho_sd",
                    label = h5("\\(\\sigma_\\rho\\)
                               Uncertainty level of \\(\\rho\\)"),
                    value = .1,
                    min = 0, max = 1
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "K",
                    label = h5("\\(K\\) Number of cluster-level covariates"),
                    value = 0, min = 0
                  ),
                  numericInput(
                    inputId = "rsq2",
                    label = h5("\\(R_2^2\\) Variance explained by
                               cluster-level covariates"),
                    value = 0,
                    min = 0, max = 1
                  )
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  selectInput(
                    inputId = "ep_al",
                    label = h5("Aims to achieve the desired..."),
                    choices = c("Expected Power",
                                "Assurance Level"),
                    selected = "Expected Power"
                  ),
                  sliderInput(
                    inputId = "power",
                    label = h5("Desired statistical power"),
                    value = .80,
                    min = 0, max = 1
                  ),
                  conditionalPanel(
                    condition = "input.ep_al == 'Assurance Level'",
                    sliderInput(
                      inputId = "al",
                      label = h5("Desired assurance level"),
                      value = .80,
                      min = 0, max = 1
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      selectInput(
                        inputId = "test",
                        label = h5("Test"),
                        choices = c("One-sided", "Two-sided"),
                        selected = "Two-sided"
                      )
                    ),
                    column(
                      width = 6,
                      numericInput(
                        inputId = "alpha",
                        label = h5("\\(\\alpha\\) Significance level"),
                        value = .05,
                        min = 0, max = 1
                      )
                    )
                  )
                ),
                column(
                  width = 6,
                  selectInput(
                    inputId = "detJn",
                    label = h5("Determine..."),
                    choices = c("Number of clusters (J)",
                                "Cluster size (n)"),
                    selected = "Number of clusters (J)"
                  ),
                  conditionalPanel(
                    condition = "input.detJn == 'Cluster size (n)'",
                    numericInput(
                      inputId = "J",
                      label = h5("Total number of clusters \\(J\\)"),
                      value = 50, min = 3
                    )
                  ),
                  conditionalPanel(
                    condition = "input.detJn == 'Number of clusters (J)'",
                    numericInput(
                      inputId = "n",
                      label = h5("Average cluster size \\(n\\)"),
                      value = 20, min = 1
                    )
                  ),
                  sliderInput(
                    inputId = "P",
                    label = h5("\\(P\\) Proportion of the clusters in treatment group"),
                    value = .5,
                    min = 0, max = 1
                  )
                )
              ),
            ),
            box(
              width = 12,
              style = 'margin:0px;',
              status = "primary",
              h3("Outputs"),
              fluidRow(
                column(
                  id = "output_plots",
                  state = "primary",
                  solidHeader = FALSE,
                  width = 12,
                  h4("Power Curves"),
                  plotOutput("plot1", height = "250px"),
                  h4("Selected Prior Distributions"),
                  plotOutput("plot2", height = "250px")
                )
              ),
              br(),
              h4("Descriptions"),
              column(
                id = "output_texts",
                state = "primary",
                solidHeader = FALSE,
                width = 12,
                textOutput("est"),
                br(),
                textOutput("Jn")
              )
            )
          )
        )
      )
    )
  )

  #### SERVER ####

  server <- function(input, output) {

    al <- reactive({
      if (input$ep_al == "Assurance Level") {
        input$al
      } else {
        NULL
      }
    })

    ret <- reactive({
      if (input$detJn == "Cluster size (n)") {
        Jn_crt2(d_est = input$d_est, d_sd = input$d_sd,
                rho_est = input$rho_est, rho_sd = input$rho_sd,
                rsq2 = input$rsq2, J = input$J, K = input$K, P = input$P,
                power = input$power, al = al(), plot = TRUE)
      } else if (input$detJn == "Number of clusters (J)") {
        Jn_crt2(d_est = input$d_est, d_sd = input$d_sd,
                rho_est = input$rho_est, rho_sd = input$rho_sd,
                rsq2 = input$rsq2, n = input$n, K = input$K, P = input$P,
                power = input$power, al = al(), plot = TRUE)
      }
    })

    output$plot1 <- renderPlot({
      do.call(gridExtra::grid.arrange, c(ret()[[1]], ncol = 2))
    })


    output$plot2 <- renderPlot({
      plot_prior(d_est = input$d_est, d_sd = input$d_sd,
                 rho_est = input$rho_est, rho_sd = input$rho_sd)
    })

    output$est <- renderText({
      paste0(
        "You provided an effect size value of ", input$d_est,
        " with an uncertainty level of ", input$d_sd,
        " and an intraclass correlation value of ", input$rho_est,
        " with a uncertainty level of ", input$rho_sd, ". "
      )
    })

    output$Jn <- renderText({
      if (input$ep_al == "Expected Power") {
        paste0("For the given inputs, a study requires J = ", ret()[[2]][1],
               " and n = ", ret()[[2]][2], " to achieve ",
               input$power*100, "% expected power. ",
               "In other words, a study with this design achieves ",
               input$power*100, "% power on average over the specified uncertainty.")
      } else {
        paste0("For the given inputs, a study requires J = ", ret()[[2]][1],
               " and n = ", ret()[[2]][2], " to achieve ",
               input$al*100, "% assurance level with ", input$power*100, "% power. ",
               "In other words, there is ", input$al*100, "% chance that a study ",
               "with this design achieves the desired statistical power. ")
      }
    })
  }

  shinyApp(ui, server)
}

