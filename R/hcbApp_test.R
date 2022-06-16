#' An Interactive Shiny App for HCB Power Analysis
#'
#' \code{hcbApp()} executes the interactive Shiny App for HCB power analysis.
#'
#' @return An interactive Shiny App to perform HCB power analysis.
#' @import shiny
#' @import shinydashboard
#' @export
#'
hcbApp_test <- function() {

  #### helper functions ####

  render_Jn_plots <- function(res) {
    do.call(gridExtra::grid.arrange, c(res$Jn_plots, ncol = 2))
  }

  render_prior_plots <- function(res) {
    do.call(gridExtra::grid.arrange, c(res$prior_plots, ncol = 2))
  }

  render_Jn <- function(res, ep_al, power, al) {
    if (ep_al == "Expected Power") {
      paste0("For the given inputs, a study requires J = ", res$Jn[1],
             " and n = ", res$Jn[2], " to achieve ",
             power*100, "% expected power. ",
             "In other words, a study with this design achieves ",
             power*100, "% power on average over the specified uncertainty.")
    } else {
      paste0("For the given inputs, a study requires J = ", res$Jn[1],
             " and n = ", res$Jn[2], " to achieve ",
             al*100, "% assurance level with ",
             power*100, "% power. ",
             "In other words, there is ", al*100, "% chance that a study ",
             "with this design achieves the desired statistical power. ")
    }
  }

  #### UI ####

  ui <- dashboardPage(
    dashboardHeader(title = "Sample Size Planning with the HCB Approach",
                    titleWidth = 500),
    dashboardSidebar(width = 150,
                     sidebarMenu(
                       menuItem(
                         "Home",
                         tabName = "home",
                         icon = icon("home")
                       ),
                       menuItem(
                         "Two-Level CRT",
                         tabName = "crt2",
                         icon = icon("users")
                       ),
                       menuItem(
                         "Two-Level MSRT",
                         tabName = "msrt2",
                         icon = icon("users")
                       )
                     )),
    dashboardBody(
      tabItems(
        tabItem(
          tabName = "home",
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
        #### Two-level CRT ####
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
                    inputId = "d_est_crt2",
                    label = h5("\\(\\delta\\) Effect size value"),
                    value = .5
                  ),
                  numericInput(
                    inputId = "rho_est_crt2",
                    label = h5("\\(\\rho\\) ICC value"),
                    value = .2,
                    min = 0, max = 1
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "d_sd_crt2",
                    label = h5("\\(\\sigma_\\delta\\)
                               Uncertatinty level of \\(\\delta\\)"),
                    value = .1, min = 0
                  ),
                  numericInput(
                    inputId = "rho_sd_crt2",
                    label = h5("\\(\\sigma_\\rho\\)
                               Uncertainty level of \\(\\rho\\)"),
                    value = .1,
                    min = 0, max = 1
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "K_crt2",
                    label = h5("\\(K\\) Number of cluster-level covariates"),
                    value = 0, min = 0
                  ),
                  numericInput(
                    inputId = "rsq2_crt2",
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
                    inputId = "ep_al_crt2",
                    label = h5("Aims to achieve the desired..."),
                    choices = c("Expected Power",
                                "Assurance Level"),
                    selected = "Expected Power"
                  ),
                  sliderInput(
                    inputId = "power_crt2",
                    label = h5("Desired statistical power"),
                    value = .80,
                    min = 0, max = 1
                  ),
                  conditionalPanel(
                    condition = "input.ep_al_crt2 == 'Assurance Level'",
                    sliderInput(
                      inputId = "al_crt2",
                      label = h5("Desired assurance level"),
                      value = .80,
                      min = 0, max = 1
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      selectInput(
                        inputId = "test_crt2",
                        label = h5("Test"),
                        choices = c("One-sided", "Two-sided"),
                        selected = "Two-sided"
                      )
                    ),
                    column(
                      width = 6,
                      numericInput(
                        inputId = "alpha_crt2",
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
                    inputId = "detJn_crt2",
                    label = h5("Determine..."),
                    choices = c("Number of clusters (J)",
                                "Cluster size (n)"),
                    selected = "Number of clusters (J)"
                  ),
                  conditionalPanel(
                    condition = "input.detJn_crt2 == 'Cluster size (n)'",
                    numericInput(
                      inputId = "J_crt2",
                      label = h5("Total number of clusters \\(J\\)"),
                      value = 50, min = 3
                    )
                  ),
                  conditionalPanel(
                    condition = "input.detJn_crt2 == 'Number of clusters (J)'",
                    numericInput(
                      inputId = "n_crt2",
                      label = h5("Average cluster size \\(n\\)"),
                      value = 20, min = 1
                    )
                  ),
                  sliderInput(
                    inputId = "P_crt2",
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
                  id = "out_plots_crt2",
                  state = "primary",
                  solidHeader = FALSE,
                  width = 12,
                  h4("Power Curves"),
                  plotOutput("Jn_plots_crt2", height = "250px"),
                  h4("Selected Prior Distributions"),
                  plotOutput("prior_plots_crt2", height = "250px")
                )
              ),
              br(),
              h4("Descriptions"),
              column(
                id = "out_texts_crt2",
                state = "primary",
                solidHeader = FALSE,
                width = 12,
                textOutput("est_crt2"),
                br(),
                textOutput("Jn_crt2")
              )
            )
          )
        ),
        #### Two-level MSRT ####
        tabItem(
          tabName = "msrt2",
          titlePanel("Design a Two-Level Multisite Randomized Trial"),
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
                    inputId = "d_est_msrt2",
                    label = h5("\\(\\delta\\) Effect size"),
                    value = .5
                  ),
                  numericInput(
                    inputId = "rho_est_msrt2",
                    label = h5("\\(\\rho\\) Intraclass correlation"),
                    value = .2,
                    min = 0, max = 1
                  ),
                  numericInput(
                    inputId = "omega_est_msrt2",
                    label = h5("\\(\\omega\\) Effect size heterogeneity"),
                    value = .2,
                    min = 0, max = 1
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "d_sd_msrt2",
                    label = h5("\\(\\sigma_\\delta\\)
                               Uncertatinty level of \\(\\delta\\)"),
                    value = .1, min = 0
                  ),
                  numericInput(
                    inputId = "rho_sd_msrt2",
                    label = h5("\\(\\sigma_\\rho\\)
                               Uncertainty level of \\(\\rho\\)"),
                    value = .1,
                    min = 0, max = 1
                  ),
                  numericInput(
                    inputId = "omega_sd_msrt2",
                    label = h5("\\(\\sigma_\\omega\\)
                               Uncertainty level of \\(\\omega\\)"),
                    value = 0,
                    min = 0, max = 1
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "K_msrt2",
                    label = h5("\\(K\\) Number of cluster-level covariates"),
                    value = 0, min = 0
                  ),
                  numericInput(
                    inputId = "rsq2_msrt2",
                    label = h5("\\(R_2^2\\) Variance explained by
                               cluster-level covariates"),
                    value = 0,
                    min = 0, max = 1
                  ),
                  numericInput(
                    inputId = "rsq1_msrt2",
                    label = h5("\\(R_1^2\\) Variance explained by
                               individual-level covariates"),
                    value = 0,
                    min = 0, max = 1
                  )
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  selectInput(
                    inputId = "ep_al_msrt2",
                    label = h5("Aims to achieve the desired..."),
                    choices = c("Expected Power",
                                "Assurance Level"),
                    selected = "Expected Power"
                  ),
                  sliderInput(
                    inputId = "power_msrt2",
                    label = h5("Desired statistical power"),
                    value = .80,
                    min = 0, max = 1
                  ),
                  conditionalPanel(
                    condition = "input.ep_al_msrt2 == 'Assurance Level'",
                    sliderInput(
                      inputId = "al_msrt2",
                      label = h5("Desired assurance level"),
                      value = .80,
                      min = 0, max = 1
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      selectInput(
                        inputId = "test_msrt2",
                        label = h5("Test"),
                        choices = c("One-sided", "Two-sided"),
                        selected = "Two-sided"
                      )
                    ),
                    column(
                      width = 6,
                      numericInput(
                        inputId = "alpha_msrt2",
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
                    inputId = "detJn_msrt2",
                    label = h5("Determine..."),
                    choices = c("Number of clusters (J)",
                                "Cluster size (n)"),
                    selected = "Number of clusters (J)"
                  ),
                  conditionalPanel(
                    condition = "input.detJn_msrt2 == 'Cluster size (n)'",
                    numericInput(
                      inputId = "J_msrt2",
                      label = h5("Total number of clusters \\(J\\)"),
                      value = 50, min = 3
                    )
                  ),
                  conditionalPanel(
                    condition = "input.detJn_msrt2 == 'Number of clusters (J)'",
                    numericInput(
                      inputId = "n_msrt2",
                      label = h5("Average cluster size \\(n\\)"),
                      value = 20, min = 1
                    )
                  ),
                  sliderInput(
                    inputId = "P_msrt2",
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
                  id = "out_plots_msrt2",
                  state = "primary",
                  solidHeader = FALSE,
                  width = 12,
                  h4("Power Curves"),
                  plotOutput("Jn_plots_msrt2", height = "250px"),
                  h4("Selected Prior Distributions"),
                  plotOutput("prior_plots_msrt2", height = "250px")
                )
              ),
              br(),
              h4("Descriptions"),
              column(
                id = "out_texts_msrt2",
                state = "primary",
                solidHeader = FALSE,
                width = 12,
                textOutput("est_msrt2"),
                br(),
                textOutput("Jn_msrt2")
              )
            )
          )
        )
      )
    )
  )

  #### SERVER ####

  server <- function(input, output) {

    #### Two-level CRT ####

    al <- reactive({
      if (input$ep_al_crt2 == "Assurance Level") input$al_crt2
      else NULL
    })

    det <- reactive({
      if (input$detJn_crt2 == "Cluster size (n)") {
        list(J = NULL, n = input$n_crt2)
      } else {
        list(J = input$J_crt2, n = NULL)
      }
    })

    res_crt2 <- reactive({
      Jn_crt2(d_est = input$d_est_crt2, d_sd = input$d_sd_crt2,
              rho_est = input$rho_est_crt2, rho_sd = input$rho_sd_crt2,
              rsq2 = input$rsq2_crt2, J = det()$J, n = det()$n,
              K = input$K_crt2, P = input$P_crt2,
              power = input$power_crt2, al = al(), plot = TRUE)
    })

    output$Jn_plots_crt2 <- renderPlot({
      render_Jn_plots(res_crt2())
    })
    output$prior_plots_crt2 <- renderPlot({
      render_prior_plots(res_crt2())
    })
    output$est_crt2 <- renderText({
      paste0(
        "You provided an effect size value of ", input$d_est_crt2,
        " with an uncertainty level of ", input$d_sd_crt2,
        " and an intraclass correlation value of ", input$rho_est_crt2,
        " with a uncertainty level of ", input$rho_sd_crt2, ". "
      )
    })
    output$Jn_crt2 <- renderText({
      render_Jn(res_crt2(), input$ep_al_crt2, input$power_crt2, al())
    })

    #### Two-Level MSRT ####

    al <- reactive({
      if (input$ep_al_msrt2 == "Assurance Level") input$al_msrt2
      else NULL
    })

    det <- reactive({
      if (input$detJn_msrt2 == "Cluster size (n)") {
        list(J = NULL, n = input$n_msrt2)
      } else {
        list(J = input$J_msrt2, n = NULL)
      }
    })

    res_msrt2 <- reactive({
      Jn_msrt2(d_est = input$d_est_msrt2, d_sd = input$d_sd_msrt2,
               rho_est = input$rho_est_msrt2, rho_sd = input$rho_sd_msrt2,
               omega_est = input$omega_est_msrt2, omega_sd = input$omega_sd_msrt2,
               rsq1 = input$rsq1_msrt2, rsq2 = input$rsq2_msrt2,
               J = det()$J, n = det()$n, K = input$K_msrt2, P = input$P_msrt2,
               power = input$power_msrt2, al = al(), plot = TRUE)
    })

    output$Jn_plots_msrt2 <- renderPlot({
      render_Jn_plots(res_msrt2())
    })
    output$prior_plots_msrt2 <- renderPlot({
      render_prior_plots(res_msrt2())
    })
    output$est_msrt2 <- renderText({
      paste0(
        "You provided an effect size of ", input$d_est_msrt2,
        " with an uncertainty level of ", input$d_sd_msrt2,
        ", an intraclass correlation of ", input$rho_est_msrt2,
        " with an uncertainty level of ", input$rho_sd_msrt2,
        ", and an effect size heteroegentiy of ", input$omega_est_msrt2,
        " with an uncertainty level of", input$omega_sd_msrt, ". "
      )
    })
    output$Jn_msrt2 <- renderText({
      render_Jn(res_msrt2(), input$ep_al_msrt2, input$power_msrt2, al())
    })
  }

  shinyApp(ui, server)
}

