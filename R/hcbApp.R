#' An Interactive Shiny App for HCB Power Analysis
#'
#' \code{hcbApp()} executes the interactive Shiny App for HCB power analysis.
#'
#' @return An interactive Shiny App to perform HCB power analysis.
#' @import shiny
#' @importFrom magrittr %>%
#' @export
#'
# hcbApp <- function() {
#   #### UI ####
#
#   ui <- fluidPage(
#     navbarPage("Hybrid Classical-Bayesian Power Analysis",
#                navbarMenu("Study Design",
#                           tabPanel("2-Level Cluster Randomized",
#                                    sidebarLayout(
#                                      sidebarPanel(
#                                        tags$style(".well {background-color: #EEFEF2;}"),
#                                        h3("Required Inputs"),
#                                        numericInput("d_est", h5("Effect size estimate"),
#                                                     value = .35),
#                                        numericInput("d_sd", h5("Uncertatinty level of the effect size estimate"),
#                                                     value = 0, min = 0),
#                                        numericInput("rho_est", h5("ICC estimate"),
#                                                     value = .15, min = 0, max = 1),
#                                        numericInput("rho_sd", h5("Uncertainty level of the ICC estimate"),
#                                                     value = 0, min = 0, max = 1),
#                                        # numericInput("r2_est", h5("Cluster-Level Covariate R2 Estimate"),
#                                        #              value = 0, min = 0, max = 1),
#                                        # numericInput("r2_sd", h5("SE of Cluster-Level Covariate R2 Estimate"),
#                                        #              value = 0, min = 0, max = 1),
#                                        # numericInput("K", h5("Number of Cluster-Level Covariate"),
#                                        #              value = 0, min = 0),
#                                        sliderInput("power", h5("Desired power level"),
#                                                    value = .80, min = 0, max = 1),
#                                        selectInput("fixedJn", "Is J or n fixed?",
#                                                    choices = c("Fixed number of clusters (J)",
#                                                                "Fixed cluster size (n)"),
#                                                    selected = "Fixed cluster size (n)"),
#                                        conditionalPanel(
#                                          condition = "input.fixedJn == 'Fixed number of clusters (J)'",
#                                          numericInput("J", h5("J"),
#                                                       value = 50, min = 3)
#                                        ),
#                                        conditionalPanel(
#                                          condition = "input.fixedJn == 'Fixed cluster size (n)'",
#                                          numericInput("n", h5("n"),
#                                                       value = 20, min = 3)
#                                        )
#                                      ),
#                                      mainPanel(
#                                        plotOutput("plot1") %>%
#                                          shinycssloaders::withSpinner(type = 6, size = .8, color = "#1DD6A2"),
#                                        br(),
#                                        br(),
#                                        textOutput("est"),
#                                        br(),
#                                        textOutput("Jn")
#                                      )
#                                    )
#                           )
#                )
#     )
#   )
#
#   #### SERVER ####
#
#   server <- function(input, output) {
#
#     output$plot1 <- output$plot2 <- output$plot3 <- renderPlot({
#       if (input$fixedJn == "Fixed number of clusters (J)") {
#         ret <- crtJn(d_est = input$d_est, d_sd = input$d_sd,
#                      rho_est = input$rho_est, rho_sd = input$rho_sd,
#                      r2_est = 0, r2_sd = 0,
#                      J = input$J,K = 0,
#                      power = input$power, plot = TRUE)
#         do.call(gridExtra::grid.arrange, c(ret[-length(ret)], ncol = 2))
#       } else if (input$fixedJn == "Fixed cluster size (n)") {
#         ret <- crtJn(d_est = input$d_est, d_sd = input$d_sd,
#                      rho_est = input$rho_est, rho_sd = input$rho_sd,
#                      r2_est = 0, r2_sd = 0,
#                      n = input$n, K = 0,
#                      power = input$power, plot = TRUE)
#         do.call(gridExtra::grid.arrange, c(ret[-length(ret)], ncol = 2))
#       }
#     })
#
#     output$est <- renderText({
#       paste0(
#         "You provided an effect size estimate of ", input$d_est,
#         " with an uncertainty level of ", input$d_sd,
#         ", an intraclass correlation estaimte of ", input$rho_est,
#         " with a uncertainty level of ", input$rho_sd, ". "
#       )
#     })
#
#     output$Jn <- renderText({
#       if (input$fixedJn == "Fixed number of clusters (J)") {
#         ret <- crtJn(d_est = input$d_est, d_sd = input$d_sd,
#                      rho_est = input$rho_est, rho_sd = input$rho_sd,
#                      r2_est = 0, r2_sd = 0,
#                      J = input$J, K = 0,
#                      power = input$power, plot = TRUE)
#         paste0("For the given inputs, it requires J = ", ret$Jn[1],
#                " and n = ", ret$Jn[2], " to achieve ",
#                input$power*100, "% expected power. ")
#       } else if (input$fixedJn == "Fixed cluster size (n)") {
#         ret <- crtJn(d_est = input$d_est, d_sd = input$d_sd,
#                      rho_est = input$rho_est, rho_sd = input$rho_sd,
#                      r2_est = 0, r2_sd = 0,
#                      n = input$n, K = 0,
#                      power = input$power, plot = TRUE)
#         paste0("For the given inputs, it requires J = ", ret$Jn[1],
#                " and n = ", ret$Jn[2], " to achieve ",
#                input$power*100, "% expected power. ")
#       }
#     })
#   }
#
#   shinyApp(ui, server)
# }
#
