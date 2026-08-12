# =============================================================================
# Persephone Shiny App — Minimal Starter Version
# Agent-Based Model for Bioarchaeology
# Anderson & DeWitte, "Known Unknowns and the Osteological Paradox"
# =============================================================================
#
# GOAL OF THIS SCRIPT:
#   This is intentionally a *small* app. It exposes exactly ONE user-changeable
#   parameter (starting population size) and uses persephone's own
#   get_default_params() for everything else. Once this works and feels
#   understandable, we can add more sliders/tabs one at a time.
#
# WHAT IS SHINY?
#   A Shiny app has two main parts:
#     1. ui     — defines what the user SEES: sliders, buttons, plots, tabs.
#     2. server — defines what HAPPENS: takes user input, runs R code,
#                  produces things (like plots) for the ui to display.
#   At the bottom, shinyApp(ui, server) glues them together and launches it.
#
# =============================================================================
library(shiny)
library(ggplot2)
library(dplyr)
library(survival)
library(persephone)   # provides Simulate_Cemetery() and get_default_params()

options(shiny.sanitize.errors = FALSE)

# =============================================================================
# Helper function: assign_age_interval()
# -----------------------------------------------------------------------------
# This isn't part of persephone — it's a small convenience function we write
# ourselves to group individual ages-at-death into bioarchaeologically
# meaningful bins (e.g. "0-1", "2-5", ...) for plotting. We define it here,
# outside of ui and server, because both the Cemetery and Survival tabs may
# want to use it, and it doesn't depend on any user input.
# =============================================================================

assign_age_interval <- function(age) {
  factor(
    case_when(
      age < 2  ~ "0-1",
      age < 6  ~ "2-5",
      age < 10 ~ "6-9",
      age < 15 ~ "10-14",
      age < 20 ~ "15-19",
      age < 30 ~ "20-29",
      age < 40 ~ "30-39",
      age < 50 ~ "40-49",
      age < 60 ~ "50-59",
      TRUE     ~ "60+"
    ),
    # Setting 'levels' explicitly keeps the bins in chronological order on
    # plots, instead of R's default alphabetical sort (which would put
    # "10-14" before "2-5").
    levels = c("0-1","2-5","6-9","10-14","15-19",
               "20-29","30-39","40-49","50-59","60+")
  )
}

# Siler survivorship function S(x): probability of surviving from birth to age x
siler_survival <- function(age, params) {
  with(params, exp(-(
    (a1 / b1) * (1 - exp(-b1 * age)) +
      a2 * age +
      (a3 / b3) * (exp(b3 * age) - 1)
  )))
}

# Life expectancy at birth: area under S(x) from 0 to upper (a large plausible max age)
life_expectancy_at_birth <- function(params, upper = 110) {
  integrate(function(x) siler_survival(x, params), lower = 0, upper = upper)$value
}

# Remaining life expectancy at a given age (e.g. age 15), i.e. e(15)
life_expectancy_at_age <- function(params, age, upper = 110) {
  s_age <- siler_survival(age, params)
  integrate(function(x) siler_survival(x, params) / s_age, lower = age, upper = upper)$value
}


# A couple of colors we'll reuse across plots, so they stay consistent.
col_no_lesion <- "#5B8DB8"   # blue
col_lesion    <- "#CC6677"   # red
col_dark      <- "#2c3e50"   # near-black
col_insample  <- "#44AA99"


# options for mortality regimes, from persephone's pre-loaded library of mortality regimes
regime_choices <- persephone::list_mortality_regimes()
# options for post-mortem preservation/taphonomic loss regimes, from persephone's pre-loaded library of taphonomy regimes
taphonomic_loss_choices <- persephone::list_taphonomy_regimes()
# =============================================================================
# UI: what the user sees
# =============================================================================
# fluidPage() is the outermost container for a simple, single-page app.
# sidebarLayout() splits the page into a narrow sidebar (controls) and a
# wider main panel (output).

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
    h4.param-header {
      color: #495057; border-bottom: 1px solid #dee2e6;
      padding-bottom: 5px; margin-top: 18px; margin-bottom: 10px; font-size: 20px;
    }
    .help-text  { font-size: 11px; color: #888; margin-top: -5px; margin-bottom: 8px; }
    .note-text  { font-size: 11px; color: #5a6a72; background: #f0f4f5;
                  border-left: 3px solid #8fb4bf; padding: 6px 8px;
                  margin-bottom: 8px; border-radius: 2px; }
    .run-btn    { margin-top: 15px; width: 100%; font-size: 16px; padding: 10px; }
    .interp     { padding: 12px; border-radius: 5px; margin-top: 10px;
                  font-size: 13px; line-height: 1.5; }
    .interp-warn  { background: #f8d7da; border: 1px solid #f5c6cb; }
    .interp-ok    { background: #d4edda; border: 1px solid #c3e6cb; }
    .summary-text { font-size: 14px; }
    .greyed-msg   { color: #888; font-style: italic; padding: 20px; }
    .conditional-block { background: #f8f9fa; border-left: 3px solid #6c9fc6;
                         padding: 8px 10px; margin-bottom: 6px; border-radius: 2px; }
  "))),

  titlePanel(div(
    "Persephone",
    tags$small(style = "color: #6c757d; display: block; font-size: 14px; margin-top: 2px;",
               "An agent-based model for bioarchaeology")
  ),
  windowTitle = "Persephone ABM"
  ),

  sidebarLayout(

    # ---------------------------------------------------------------------
    # SIDEBAR: user controls go here.
    # ---------------------------------------------------------------------
    sidebarPanel(
      width = 4,
      # ---- Run Simulation -----------------------------------------------
      # We use this (rather than running the simulation automatically every
      # time the slider moves) so the simulation only re-runs when the user
      # is ready — agent-based models can take a few seconds to run, so we
      # don't want to re-run on every tiny slider nudge.
      actionButton(
        inputId = "run",
        label   = "Run Simulation",
        class   = "btn-primary",   # just a Bootstrap CSS class for blue styling
      ),
      
      ### Max Years to Run Simulation ###
      sliderInput("max_years", "Simulation duration (years)",
                  value = 100, min = 10, max = 1000, step = 10),
      div(class = "help-text", "If the population does not crash, the simulation will run until you tell it to stop."),
      # ---- Population -------------------------------------------------------
      h4("Population", class = "param-header"),
      
      ### A slider for starting population size. ###
      # sliderInput(inputId, label, min, max, value, step)
      #   - inputId = "pop0_size" is the name we'll use in server() to read
      #     whatever value the user has selected, via input$pop0_size.
      sliderInput(
        inputId = "pop0_size",
        label   = "Starting population size",
        min     = 50,
        max     = 5000,
        value   = 1000,   # matches persephone's default
        step    = 50
      ),
      
      ### An option for specifying a stationary or non-stationary (but stable) population. ###
      h4("", class = "param-header"),
      div(class = "help-text", "To model a growing or shrinking population, un-check the box below."),
      checkboxInput(inputId = "nonstationarity",
                    label = "Stationary Population (pop. growth = 0)",
                    value = TRUE),

      conditionalPanel(
        condition = "input.nonstationarity == false",
        div(class = "help-text", "To model a non-stationary population, you need to specify how the population is growing."),
        
        ### Total Fertility Rate ###
        h4("Fertility", class = "param-header"),
        sliderInput("tfr", "Total fertility rate (TFR)", value = 2, min = 0, max = 12, step = 0.1),
        div(class = "help-text", "Average number of children born to a woman in this population.")
        ),
      
      # ---- Mortality --------------------------------------------------------
      h4("Mortality", class = "param-header"),
      selectInput(inputId = "mortality_regime", 
                  label = "Mortality regime",
                  choices = names(regime_choices), 
                  selected = "CoaleDemenyWestF5"),
      
      # ---- Skeletal Lesions -------------------------------------------------
      h4("Morbidity", class = "param-header"),
      checkboxInput("model_lesions", "Generate skeletal lesions", value = FALSE),
      div(class = "help-text", "If you check this box, you will need to specify lesion formation dynamics"),
      conditionalPanel(
        condition = "input.model_lesions == true",
        h4("Skeletal Lesions", class = "param-header"),
        fluidRow(
          column(6, numericInput("window_opens",  "Window opens",  0, min = 0, max = 80)),
          column(6, numericInput("window_closes", "Window closes", 9, min = 0, max = 80))
        ),
        div(class = "help-text", "Age range in which skeletal lesions can form."),
            sliderInput("lesion_formation_rate", "Annual probability of encountering lesion-causing conditions",
                        min = 0.01, 
                        max = 1, 
                        value = 0.05, # default value
                        step = 0.01)
      ),
      
      # ---- Post-mortem processes --------------------------------------------
      h4("Post-mortem processes", class = "param-header"),
      checkboxInput("model_postmortem", "Simulate post-mortem processes", value = FALSE),
      div(class = "help-text", "Check this box to add burial practices, preservation loss, and observation error during data collection to your simulated bioarchaeological data set."),
      conditionalPanel(
        condition = "input.model_postmortem == true",
        h4("Post-Mortem Processes", class = "param-header"),
      numericInput("deposition_param", "Minimum deposition age (years)",
                   0, min = 0, max = 10, step = 1),
      div(class = "help-text",
          "Individuals below this age are excluded from the cemetery."),
      
      checkboxInput("enable_taphonomy", "Apply taphonomic loss", FALSE),
      conditionalPanel(
        condition = "input.enable_taphonomy == true",
        div(class = "conditional-block",
            selectInput("loss_strength", "Preservation Loss",
                        choices = c("None" = "no_decay",
                                    "Mild" = "weak_decay", 
                                    "Moderate" = "moderate_decay",
                                    "Severe" = "strong_decay"),
                        selected = "moderate"))
      ),
      
      checkboxInput("age_noise", "Add skeletal age estimation error", FALSE),
      div(class = "help-text",
          "Adds stochastic noise (and bias, at older ages) to estimated age-at-death to simulate osteological assessment error.")
      ),
      
      # ---- Make Replicable (set seed) --------------------------------------------------------
      h4("Simulation", class = "param-header"),
      # Set seed for reproducible results:
      numericInput("seed", "Random seed (blank = random)", 1, min = 1),
    ),
      

    # ---------------------------------------------------------------------
    # MAIN PANEL: output goes here, organized into tabs.
    # tabsetPanel() creates a set of tabs; each tabPanel() is one tab.
    # ---------------------------------------------------------------------
    mainPanel(
      # ---- Tab 1: Parameter Explanations -----------------------------
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel(
          "Model Parameters",
        br(),
        uiOutput("param_explanation")
        ),

        # ---- Tab 2: Population Dynamics ------------------------------------------
          tabPanel("Population Dynamics",
                   br(),
                   uiOutput("pop_dynamics_ui")
          ),
        
        # ---- Tab 3: Cemetery -----------------------------------
        tabPanel(
          "Cemetery",
          br(),  
          uiOutput("cemetery_ui")
        ),
        # ---- Tab 4: Survival Analysis -----------------------------------
        tabPanel(
           "Survival Analysis",
           br(),
           uiOutput("survival_ui")
          # plotOutput("plot_km"),          # Kaplan-Meier survival curves
          # br(),
          # verbatimTextOutput("logrank_text")  # plain-text log-rank test result
        ),
        # ---- Tab 4: Multi-State Model Analysis -----------------------------------
        
        # ---- Tab 5: Survival Analysis -----------------------------------
        
        tabPanel("Data",
                 br(),
                 fluidRow(
                   column(6,
                          h4("Cemetery (Individual Outcomes)"),
                          downloadButton("dl_cemetery", "Download CSV"),
                          div(style = "max-height: 500px; overflow-y: auto; margin-top: 10px;",
                              tableOutput("table_cemetery"))
                   ),
                   column(6,
                          h4("Annual Census"),
                          downloadButton("dl_census", "Download CSV"),
                          div(style = "max-height: 500px; overflow-y: auto; margin-top: 10px;",
                              tableOutput("table_census"))
                   )
                 )
        )
      #   tabPanel("Parameters",
      #            br(),
      #            verbatimTextOutput("sim_params"))
       )
    )
  )
)
    



# =============================================================================
# SERVER: what happens when the user interacts with the app
# =============================================================================
# server is a function that takes (input, output, session). We mostly use
# input (to read what the user picked) and output (to send results back to
# the ui). 'session' isn't needed for this simple app but Shiny expects the
# argument to be there.

server <- function(input, output, session) {

  # ---------------------------------------------------------------------
  # Step 1: Get persephone's default parameter list.
  # ---------------------------------------------------------------------
  # get_default_params() returns a named list where every element is one
  # argument to Simulate_Cemetery() (e.g. params$tfr, params$mortality_regime,
  # etc). We grab this ONCE, outside of any reactive code, since the defaults
  # themselves don't depend on user input — only pop0_size will be
  # overridden later, every time the user clicks "Run Simulation".
  params <- get_default_params()
  
  # to create the set of parameter values that will be used by the ABM below,
  # overwrite specified parameters with the values chosen in the control panel.
  sim_params <- reactive({
    # set a seed for start of stochastic processes. Necessary for replicable results. 
    set.seed(input$seed)
    
    # Take persephone's defaults...
    sim_params <- params
    
    # ...and overwrite just the one value we're letting the user control.
    sim_params$pop0_size <- input$pop0_size
    sim_params$max_years <- input$max_years
    
    sim_params$age_structured <- input$nonstationarity == FALSE
    # to model a dynamic population:
    if (sim_params$age_structured) {
      sim_params$tfr <- input$tfr
    }
    # to specify mortality schedule:
    sim_params$mortality_regime <- persephone::list_mortality_regimes()[[input$mortality_regime]]
    
    # to model lesions:
    if(input$model_lesions){
      sim_params$lesion_formation_rate <- input$lesion_formation_rate
      sim_params$lesion_formation_window <- c(input$window_opens, input$window_closes)
    }
    # to model post-mortem processes:
    if(input$model_postmortem){
      sim_params$deposition_param <- input$deposition_param
      sim_params$taphonomy_regime <- persephone::list_taphonomy_regimes()[input$loss_strength]
      sim_params$loss_strength <- names(sim_params$taphonomy_regime)
      sim_params$age_noise <- input$age_noise
    }
    sim_params
  })
  
  # output$sim_params <-renderUI({
  #   sim_params()
  # })
  # Calculate life expectancy stats for the chosen mortality regime
  mortality_stats <- reactive({
    regime_params <- sim_params()$mortality_regime
    
    e0                     <- life_expectancy_at_birth(regime_params)
    pct_die_before_15      <- (1 - siler_survival(15, regime_params)) * 100
    mean_age_death_after15 <- 15 + life_expectancy_at_age(regime_params, 15)
    
    list(
      e0                     = e0,
      pct_die_before_15      = pct_die_before_15,
      mean_age_death_after15 = mean_age_death_after15
    )
  })
  
  # Explain to the user what their parameter choices mean
  output$param_explanation <- renderUI({
    p <- sim_params()
    m <- mortality_stats()
    
    bold <- function(x) paste0("<b>", x, "</b>")
    italic <- function(x) paste0("<i>", x, "</i>")
    
    HTML(paste0(
      "<p><h4>Welcome to the god's eye view</h4></p>
      <p>Persephone is an agent-based model designed to generate, visualize, and analyze simulated bioarchaeological data. </p>
      <p> In the tiny world of the model, we can explore what happens to the observable data when we change a single aspect of the processes upstream from our observations. While we cannot directly observe conditions in the past, we can simulate many possible pasts and the data set that results from each possible scenario. </p>",
      "<p>To build the model, a series of decisions (outlined below) are made about how the world of the model operates. You, in choosing the parameter values in the panel to the left, are deciding the characteristics of the population and the nature of the experiences that the individuals in this population will be subjected to. </p>", 
      "<br>",
      "<p><h4>In the current model scenario: </h4></p>",
      "<p>This simulated population begins with ", bold(format(p$pop0_size, big.mark = ",")),
      " agents/individuals, modeled over ", bold(p$max_years), " years.</p>", 
      if (p$age_structured) {
        paste0("<p>This is a stable population with an average total fertility rate of ",
               bold(p$tfr), ".</p>")
      } else {
        "<p>This is a cohort model, which can be interpreted as a stationary population (zero population growth).</p>"
        },
      "<p> Before adding in any additional sources of mortality risk, this population's mortality schedule is described by a ", bold(input$mortality_regime), " mortality regime. </p>",
      "<p>This mortality curve describes an imaginary population with a life expectancy at birth (e<sub>0</sub>) of ",
      bold(round(m$e0, 1)), " years.</p>",
      
      "<p>", bold(round(m$pct_die_before_15, 1)), "% of the population does not survive to age 15.</p>",
      
      "<p>For those who do survive past age 15, mean age at death is ",
      bold(round(m$mean_age_death_after15, 1)), " years.</p>",
      
       if(input$model_lesions){
         paste0("<p>Individuals in this population have a ", bold(p$lesion_formation_rate), " annual probability of encountering skeletal-lesion-causing conditions each year between the ages of ", bold(input$window_opens), " and ", bold(input$window_closes), " years. </p>")
       },
      if(p$deposition_param > 0){
        paste0("<p>Cultural influences determine who ends up in a cemetery. In this population, individuals who die at ages younger than ", bold(p$deposition_param), " years old are not deposited in the same burial context as the rest of the population.")
      },
      if(input$enable_taphonomy && input$loss_strength != "no_decay"){
        paste0("<p>Post-burial preservation of skeletal remains is rarely perfect. In this population the extent of post-burial/taphonomic loss is ",
               bold(names(which(c("None"="no_decay","Mild"="weak_decay","Moderate"="moderate_decay","Severe"="strong_decay") == input$loss_strength))),
               ". </p>")
      },
      if(p$age_noise){
        paste0("<p>Data are the result of observation and measurement, and these are also rarely perfect. Some observation error has been added to the age-at-death data. But you will still be able to compare these simulated osteological age-at-death estimates to the true age-at-death values for individuals in this population. </p>")
      },
      "<br>",
      italic("<p><b> A Note on Mortality:</b> Currently, each mortality regime option is a Siler mortality hazard function fit to one of Coale and Demeny’s West series of model life tables for females (Coale and Demeny, 1966; Gage and Dyke, 1986). Coale and Demeny’s regional model life tables were developed based on patterns of historic European mortality and have been widely used for estimating mortality in situations with sparse data.</p>"),
      "<br>",
      italic("<p><b> A Note on Taphonomy:</b> Post-burial preservation loss is most marked in the oldest and youngest individuals in a cemetery. This pattern is qualitatively similar to mortality hazards, and has been operationalized here using Siler functions, like the mortality regimes.</p>")
    ))
  })
  # ---------------------------------------------------------------------
  # Step 2: Run the simulation when the user clicks "Run Simulation".
  # ---------------------------------------------------------------------
  # eventReactive(trigger, { code }) creates a "reactive expression" that
  # only re-runs its code when 'trigger' changes — here, input$run (the
  # button). This is different from a normal reactive(), which would
  # re-run any time ANY input it depends on changes.
  #
  # The result of this block (called by typing sim() elsewhere in server)
  # is the list that Simulate_Cemetery() returns: a list with
  # $individual_outcomes (one row per simulated skeleton) and
  # $annual_census (population size each year).
  sim <- eventReactive(input$run, {
    validate(
      need(input$pop0_size >= 50, "Population size must be at least 50."),
      need(
        !input$model_lesions || input$window_closes >= input$window_opens,
        "Formation window: 'closes' must be >= 'opens'."
      )
    )
    # do.call() lets us call Simulate_Cemetery() using a *list* of named
    # arguments, rather than typing every argument name by hand. Since
    # sim_params already has the correct names (dx, max_years, pop0_size,
    # tfr, mortality_regime, ...) matching Simulate_Cemetery's arguments,
    # this is equivalent to writing out all ~20 arguments explicitly.
    #
    # withProgress() just shows a small progress bar/message in the app
    # while the simulation runs, so the user knows something is happening.
    withProgress(message = "Running simulation...", value = 0.3, {
      print(str(sim_params()))
      result <- do.call(Simulate_Cemetery, sim_params())
      setProgress(1)
    })

    result
  })

  # Indicate whether simulation models lesions
  has_lesions <- reactive({
    req(sim())
    d <- sim()$individual_outcomes
    ("lesion" %in% names(d) || "Lesion" %in% names(d))
  })
  
  # Indicate whether simulation models post-mortem processes
  has_postmortem <- reactive({
    req(sim())
    d <- sim()$individual_outcomes %>%
      filter(!is.na(in_sample))
    sum(d$in_sample) < nrow(d) || "estimated_age" %in% names(d)
  })
  
  # ---------------------------------------------------------------------
  # Step 3: Build a "tidy" version of the cemetery data for plotting.
  # ---------------------------------------------------------------------
  # This is a normal reactive() (not eventReactive()), so it automatically
  # re-runs whenever sim() changes (i.e. right after a new simulation run).
  # We do small cleanup here once, instead of repeating it in every plot.
  cemetery <- reactive({

    # req(sim()) tells Shiny "don't run any of this until sim() has a
    # value" — i.e. until the user has clicked "Run Simulation" at least
    # once. Without this, Shiny would try to draw empty plots on startup
    # and throw confusing errors.
    req(sim())

    d <- sim()$individual_outcomes

    # persephone's default lesion-formation mode is 'annual_exposure'
    # (see get_default_params() output), so this simulation WILL include
    # a lesion column. The column may be lowercase ('lesion') or
    # capitalized ('Lesion') depending on the package version, so we
    # standardize it here to avoid that tripping up the plotting code.
    if ("lesion" %in% names(d) && !"Lesion" %in% names(d)) {
      names(d)[names(d) == "lesion"] <- "Lesion"
    }
    if ("age" %in% names(d) && !"Age" %in% names(d)) {
      names(d)[names(d) == "age"] <- "Age"
    }

    # Add the age-interval bins and a readable Yes/No lesion label,
    # using the helper function we defined above the UI/server code.
    d$Age_Interval  <- assign_age_interval(d$Age)
    if('estimated_age' %in% names(d)){
      d$estimated_age <- as.numeric(d$estimated_age)
      d$Age_Interval_Est <- assign_age_interval(d$estimated_age)
    }
    if('Lesion' %in% names(d)){
    d$Lesion_Status <- 
      factor(
      ifelse(d$Lesion == 1, "Present", "Absent"),
      levels = c("Absent", "Present")
      ) 
    } else{
      d$Lesion_Status <- NA
      }

    d
  })

  
    # ===========================================================================
    # Tab 2: Population Dynamics
    # ===========================================================================
  output$pop_dynamics_ui <- renderUI({
    if (is.null(sim())) return(div(class = "greyed-msg",
                                   "Click 'Run Simulation' to see population dynamics."))
    if (isTRUE(sim_params()$age_structured)) {
      tagList(
        fluidRow(column(12, plotOutput("plot_pop_size", height = "420px"))),
        div(class = "help-text", style = "padding: 10px;",
            "Annual census of the living population across simulation time.")
      )
    } else {
      tagList(fluidRow(
        column(6, plotOutput("plot_alive",      height = "400px")),
        column(6, plotOutput("plot_lesion_pct", height = "400px"))
      ))
    }
  })   
  

    output$plot_pop_size <- renderPlot({
      req(sim())
      census <- as.data.frame(sim()$annual_census)
      ggplot(census, aes(x = Time, y = n)) +
        geom_area(fill = col_no_lesion, alpha = 0.3) +
        geom_line(color = col_dark, linewidth = 1) +
        labs(title = "Living Population Size Over Time",
             x = "Simulation Year", y = "Number Alive") +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    })

    output$plot_alive <- renderPlot({
      req(sim())
      census <- as.data.frame(sim()$annual_census)
      ggplot(census, aes(x = Age, y = Alive)) +
        geom_area(fill = col_no_lesion, alpha = 0.3) +
        geom_line(color = col_dark, linewidth = 1) +
        labs(title = "Surviving Cohort Over Time",
             x = "Age (years)", y = "Number Alive") +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    })

    output$plot_lesion_pct <- renderPlot({
      req(sim())
      if (!has_lesions()) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "No lesions modeled in this simulation.",
                   color = "#888", size = 5, hjust = 0.5) + theme_void())
      }
      census <- as.data.frame(sim()$annual_census)
      validate(need("Lesion_perc" %in% names(census), "No lesion data available."))
      surv <- census %>% filter(!is.na(Lesion_perc))
      ggplot(surv, aes(x = Age, y = Lesion_perc)) +
        geom_area(fill = col_lesion, alpha = 0.2) +
        geom_line(color = col_lesion, linewidth = 1) +
        labs(title = "Lesion Prevalence Among Living (Cohort)",
             x = "Age (years)", y = "% with Lesions") +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        coord_cartesian(ylim = c(0, NA))
    })

  # ---------------------------------------------------------------------
  # Tab 3: Cemetery Age-at-Death Distribution and Lesion Distribution
  # ---------------------------------------------------------------------
    output$cemetery_ui <- renderUI({
      validate(need(sim(), "Click 'Run Simulation' to begin."))
      
      show_insample <- has_postmortem()   # already means "in_sample differs from full sample"
      show_lesion   <- has_lesions()
      
      age_row <- if (show_insample) {
        fluidRow(
          column(6, plotOutput("plot_cem_age_full",     height = "400px")),
          column(6, plotOutput("plot_cem_age_insample",  height = "400px"))
        )
      } else {
        fluidRow(column(12, plotOutput("plot_cem_age_full", height = "400px")))
      }
      
      lesion_row <- if (show_lesion) {
        if (show_insample) {
          fluidRow(
            column(6, plotOutput("plot_cem_lesion_full",    height = "350px")),
            column(6, plotOutput("plot_cem_lesion_insample", height = "350px"))
          )
        } else {
          fluidRow(column(12, plotOutput("plot_cem_lesion_full", height = "350px")))
        }
      } else NULL
      
      tagList(age_row, lesion_row,
     wellPanel(h4("Mean ages at death", class = "param-header"),
                               verbatimTextOutput("cem_summary"))
      )
    })
    
    
    cemetery_display <- reactive({
      req(sim())
      d <- cemetery()
      list(
        data     = d,
        has_est  = isTRUE(sim_params()$age_noise) && "Age_Interval_Est" %in% names(d)
      )
    })
    
    plot_age_bar <- function(data, interval_col, fill_lesion = FALSE,
                             x_lab = "Age at Death", title = "Age-at-Death Distribution") {
      p <- ggplot(data, aes(x = .data[[interval_col]]))
      if (fill_lesion) {
        p <- p + aes(fill = Lesion_Status) +
          geom_bar(position = "stack", width = 0.8) +
          scale_fill_manual(values = c("Absent" = col_no_lesion, "Present" = col_lesion), name = "Lesion")
      } else {
        p <- p + geom_bar(fill = col_no_lesion, width = 0.8)
      }
      p + labs(title = title, x = x_lab, y = "Count") +
        theme_bw(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "bottom")
    }
    
    plot_lesion_pct_bar <- function(data, interval_col, x_lab = "Age at Death",
                                    title = "Lesion Prevalence by Age Group") {
      prev <- data %>%
        group_by(.data[[interval_col]]) %>%
        summarise(Pct = sum(Lesion == 1, na.rm = TRUE) / n() * 100, n = n(), .groups = "drop")
      ggplot(prev, aes(x = .data[[interval_col]], y = Pct)) +
        geom_col(fill = col_lesion, alpha = 0.85, width = 0.8) +
        geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3.2, color = "#555") +
        labs(title = title, x = x_lab, y = "% with Lesions") +
        theme_bw(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5, face = "bold")) +
        coord_cartesian(ylim = c(0, max(prev$Pct, na.rm = TRUE) * 1.2))
    }
    
    
    output$plot_cem_age_full <- renderPlot({
      disp <- cemetery_display()
      # Full sample always uses TRUE age, regardless of age_noise
      plot_age_bar(disp$data, "Age_Interval", fill_lesion = has_lesions(),
                   x_lab = "Age at Death", title = "Age-at-Death Distribution (Full Sample)")
    })
    
    output$plot_cem_age_insample <- renderPlot({
      disp <- cemetery_display()
      d <- disp$data %>% filter(in_sample)
      # Recovered sample uses ESTIMATED age, if age_noise is on
      interval_col <- if (disp$has_est) "Age_Interval_Est" else "Age_Interval"
      x_lab        <- if (disp$has_est) "Estimated Age at Death" else "Age at Death"
      plot_age_bar(d, interval_col, fill_lesion = has_lesions(),
                   x_lab = x_lab, title = "Age-at-Death Distribution (Recovered Sample)")
    })
    
    output$plot_cem_lesion_full <- renderPlot({
      req(has_lesions())
      disp <- cemetery_display()
      plot_lesion_pct_bar(disp$data, "Age_Interval", "Age at Death",
                          "Lesion Prevalence (Full Sample)")
    })
    
    output$plot_cem_lesion_insample <- renderPlot({
      req(has_lesions())
      disp <- cemetery_display()
      d <- disp$data %>% filter(in_sample)
      interval_col <- if (disp$has_est) "Age_Interval_Est" else "Age_Interval"
      x_lab        <- if (disp$has_est) "Estimated Age at Death" else "Age at Death"
      plot_lesion_pct_bar(d, interval_col, x_lab, "Lesion Prevalence (Recovered Sample)")
    })
  
  
  # Tab 3, Plot 3: Cemetery summary
  output$cem_summary <- renderPrint({
    validate(need(sim(), ""))
    d <- cemetery()
  cat(sprintf("Total individuals in simulation: %d\n", nrow(d)))
  if (has_postmortem()) {
    n_in <- sum(d$in_sample, na.rm = TRUE)
    cat(sprintf("In recovered sample (post-mortem filters): %d (%.1f%%)\n",
                n_in, n_in / nrow(d) * 100))
  }
  cat(sprintf("Mean age at death: %.1f years\n", mean(d$Age, na.rm = TRUE)))
  if (sim_params()$age_noise){
    cat(sprintf("Observed mean age at death: %.1f years\n", mean(d$estimated_age, na.rm = TRUE)))
  }
  if (has_lesions()) {
    with_les <- sum(d$Lesion == 1, na.rm = TRUE)
    cat(sprintf("With lesions: %d (%.1f%% of cemetery)\n",
                with_les, with_les / nrow(d) * 100))
    if (with_les > 0) {
      cat(sprintf("  Mean age at death — with lesion:    %.1f years\n",
                  mean(d$Age[d$Lesion == 1], na.rm = TRUE)))
      cat(sprintf("  Mean age at death — without lesion: %.1f years\n",
                  mean(d$Age[d$Lesion == 0], na.rm = TRUE)))
      if (sim_params()$age_noise){
        cat(sprintf("Observed mean age at death — with lesion:    %.1f years\n",
                    mean(d$estimated_age[d$Lesion == 1], na.rm = TRUE)))
        cat(sprintf("Observed mean age at death — with outlesion:    %.1f years\n",
                    mean(d$estimated_age[d$Lesion == 0], na.rm = TRUE)))
      }
    }
  }
  })

  # ---------------------------------------------------------------------
  # Tab 4, Plot: Kaplan-Meier survival curves, split by lesion status
  # ---------------------------------------------------------------------
  # This compares the survival experience of individuals with vs. without
  # lesions — the core comparison underlying the osteological paradox.
       output$survival_ui <- renderUI({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    
        full_panel <- tagList(
          h4("Full Mortality Sample", style = "color: #2c3e50;"),
          div(class = "help-text", "All individuals who died in the simulation."),
          fluidRow(
            column(7, plotOutput("plot_km_full", height = "420px")),
            column(5,
              wellPanel(
                h4("Age Filter", class = "param-header"),
                sliderInput("min_age_full", "Minimum age to include", 0, 25, 0, 1),
                div(class = "help-text",
                    "Exclude deaths within the lesion formation window to isolate selective mortality effects.")
              ),
              wellPanel(h4("Log-Rank Test", class = "param-header"),
                        verbatimTextOutput("logrank_full")),
              uiOutput("interp_full")
            )
          )
        )

        if (has_postmortem()) {
          tagList(
            full_panel,
            hr(),
            h4("Recovered Sample Only", style = "color: #44AA99;"),
            div(class = "help-text",
                "Only individuals retained after deposition and taphonomy filters (in_sample == TRUE). ",
                "Compare with the full sample above to see how post-mortem processes bias inference."),
            fluidRow(
              column(7, plotOutput("plot_km_insample", height = "420px")),
              column(5,
                wellPanel(
                  h4("Age Filter", class = "param-header"),
                  sliderInput("min_age_insample", "Minimum age to include", 0, 25, 0, 1)
                ),
                wellPanel(h4("Log-Rank Test", class = "param-header"),
                          verbatimTextOutput("logrank_insample"))#,
                # uiOutput("interp_insample")
              )
            )
          )
        } else {
          full_panel
        }
      })

  # Shared KM helpers
  km_data <- function(data, min_age, age_col = "Age") {
    d <- data %>% filter(.data[[age_col]] >= min_age) %>% mutate(Dead = 1)
    if (nrow(d) < 10) return(NULL)   # only a real "too few individuals" bailout
    
    has_lesion_col    <- "Lesion" %in% names(d)
    has_lesion_groups <- has_lesion_col && length(unique(na.omit(d$Lesion))) == 2
    
    form <- if (has_lesion_groups) {
      as.formula(paste0("Surv(", age_col, ", Dead) ~ Lesion"))
    } else {
      as.formula(paste0("Surv(", age_col, ", Dead) ~ 1"))
    }
    fit <- survfit(form, data = d)
    s   <- summary(fit)
    
    df <- data.frame(
      time     = s$time,
      survival = s$surv,
      group    = if (has_lesion_groups) {
        ifelse(grepl("=0", as.character(s$strata)), "No Lesion", "Has Lesion")
      } else {
        "Overall"
      }
    )
    starts <- df %>%
      group_by(group) %>%
      summarise(time = max(min(time) - 0.5, 0), .groups = "drop") %>%
      mutate(survival = 1.0)
    dplyr::bind_rows(starts, df) %>% arrange(group, time)
  }
  
  km_ggplot <- function(sdf, min_age, title_sfx = "", colors = NULL, x_lab = "Age at Death") {
    if (is.null(sdf)) return(ggplot() +
                               annotate("text", x = 0.5, y = 0.5,
                                        label = "Too few individuals or only one lesion group after filtering.",
                                        color = "#888", size = 4.5, hjust = 0.5) + theme_void())
    if (is.null(colors)) colors <- c("No Lesion" = col_dark, "Has Lesion" = col_lesion, "Overall" = col_dark)
    ttl <- paste0("Kaplan-Meier Survival Curves", title_sfx)
    if (min_age > 0) ttl <- paste0(ttl, "  (ages \u2265 ", min_age, ")")
    ggplot(sdf, aes(x = time, y = survival, color = group)) +
      geom_step(linewidth = 1.1) +
      scale_color_manual(values = colors) +
      labs(title = ttl, x = x_lab, y = "Survival Probability", color = "") +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom", legend.text = element_text(size = 12))
  }
  
  
  logrank_res <- function(data, age_col, min_age) {
    d <- data %>% filter(.data[[age_col]] >= min_age) %>% mutate(Dead = 1)
    if (nrow(d) < 10 || length(unique(d$Lesion)) < 2) return(list(chisq = NA, p = NA))
    
    form <- as.formula(paste0("Surv(", age_col, ", Dead) ~ Lesion"))
    test <- survdiff(form, data = d)
    p    <- 1 - pchisq(test$chisq, df = length(test$n) - 1)
    list(chisq = test$chisq, p = p)
  }

      # interp_html <- function(p, rmr) {
      #   if (is.na(p)) return(NULL)
      #   sig <- p < 0.05
      #   if      (rmr == 1 && sig)  { cls <- "interp interp-warn"; msg <- paste(
      #     "FALSE POSITIVE: Significant difference detected, but lesions have no mortality effect",
      #     "(hazard multiplier = 1). This is the osteological paradox — age structure alone creates",
      #     "the signal. Try raising the minimum age above the formation window.") }
      #   else if (rmr == 1 && !sig) { cls <- "interp interp-ok";   msg <-
      #     "Correct inference: no significant difference, and lesions truly have no mortality effect." }
      #   else if (rmr > 1 && sig)   { cls <- "interp interp-ok";   msg <- paste0(
      #     "Significant difference detected — lesions do increase mortality (hazard multiplier = ",
      #     rmr, "). Check that the direction of effect in the survival curves matches the true effect.") }
      #   else                       { cls <- "interp interp-warn"; msg <- paste0(
      #     "FALSE NEGATIVE: No significant difference despite a true mortality effect (hazard multiplier = ",
      #     rmr, "). Selective mortality or age confounding may be masking the signal.") }
      #   div(class = cls, HTML(msg))
      # }
      
      
      output$plot_km_full <- renderPlot({
        req(sim())
        min_a <- if (!is.null(input$min_age_full)) input$min_age_full else 0
        # Full sample always uses TRUE age, regardless of age_noise
        km_ggplot(km_data(cemetery(), min_a, age_col = "Age"), min_a,
                  x_lab = "Age at Death")
      })
      
      output$plot_km_insample <- renderPlot({
        req(sim(), has_postmortem())
        d     <- cemetery() %>% filter(in_sample == TRUE)
        min_a <- if (!is.null(input$min_age_insample)) input$min_age_insample else 0
        
        use_est <- isTRUE(sim_params()$age_noise) && "estimated_age" %in% names(d)
        age_col <- if (use_est) "estimated_age" else "Age"
        x_lab   <- if (use_est) "Estimated Age at Death" else "Age at Death"
        
        km_ggplot(km_data(d, min_a, age_col = age_col), min_a,
                  title_sfx = "",
                  colors = c("No Lesion" = col_dark, "Has Lesion" = col_insample, "Overall" = col_dark),
                  x_lab = x_lab)
      })
      
      output$logrank_full <- renderPrint({
        req(sim(), has_lesions())
        d     <- cemetery() %>% filter(in_sample == TRUE)
        min_a <- if (!is.null(input$min_age_full)) input$min_age_full else 0
        age_col <-  if('estimated_age' %in% names(d)) "estimated_age" else "Age"
        res   <- logrank_res(cemetery(), age_col, min_a)
        if (is.na(res$p)) { cat("Too few individuals or one group only."); return() }
        cat(sprintf("Chi-squared: %.2f\n", res$chisq))
        cat(sprintf("p-value:     %.4f\n", res$p))
        cat(sprintf("\n%s at alpha = 0.05", ifelse(res$p < 0.05, "SIGNIFICANT", "Not significant")))
      })
      
      # output$interp_full <- renderUI({
      #   req(sim(), has_lesions())
      #   min_a <- if (!is.null(input$min_age_full)) input$min_age_full else 0
      #   interp_html(logrank_res(cemetery(), min_a)$p, input$lesion_related_hazard)
      # })
      # 
      
      output$logrank_insample <- renderPrint({
        req(sim(), has_lesions(), has_postmortem())
        d     <- cemetery() %>% filter(in_sample == TRUE)
        min_a <- if (!is.null(input$min_age_insample)) input$min_age_insample else 0
        age_col <-  if('estimated_age' %in% names(d)) "estimated_age" else "Age"
        res   <- logrank_res(d, age_col, min_a)
        if (is.na(res$p)) { cat("Too few individuals or one group only."); return() }
        cat(sprintf("Chi-squared: %.2f\n", res$chisq))
        cat(sprintf("p-value:     %.4f\n", res$p))
        cat(sprintf("\n%s at alpha = 0.05", ifelse(res$p < 0.05, "SIGNIFICANT", "Not significant")))
      })
      # output$interp_insample <- renderUI({
      #   req(sim(), has_lesions(), has_postmortem())
      #   d     <- cemetery() %>% filter(in_sample == TRUE)
      #   min_a <- if (!is.null(input$min_age_insample)) input$min_age_insample else 0
      #   interp_html(logrank_res(d, min_a)$p, input$lesion_related_hazard)
      #  })
  # ---------------------------------------------------------------------
  # Tab 4, Text: Log-rank test (formal statistical comparison of the
  # two survival curves above)
  # ---------------------------------------------------------------------
  # renderPrint({...}) captures whatever gets printed (via cat() or print())
  # inside the block and displays it as plain text in the app — this pairs
  # with verbatimTextOutput("logrank_text") in the ui.
  output$logrank_text <- renderPrint({
    req(sim())

    d <- cemetery()
    d$Dead <- 1

    validate(need(
      length(unique(d$Lesion)) == 2,
      "Need both lesion groups present to run a log-rank test."
    ))

    test <- survdiff(Surv(Age, Dead) ~ Lesion, data = d)
    p    <- 1 - pchisq(test$chisq, df = length(test$n) - 1)

    cat(sprintf("Chi-squared: %.2f\n", test$chisq))
    cat(sprintf("p-value:     %.4f\n", p))
    cat(sprintf("\n%s at alpha = 0.05",
                ifelse(p < 0.05, "SIGNIFICANT", "Not significant")))
  })


  # ===========================================================================
  # Tab 5: Data
  # ===========================================================================

  output$table_cemetery <- renderTable({
    validate(need(sim(), "Run a simulation to see data."))
    sim()$individual_outcomes
  }, digits = 1)

  output$table_census <- renderTable({
    validate(need(sim(), ""))
    as.data.frame(sim()$annual_census)
  }, digits = 1)

  output$dl_cemetery <- downloadHandler(
    filename = function() paste0("persephone_cemetery_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) { req(sim()); write.csv(sim()$individual_outcomes, file, row.names = FALSE) }
  )
}

shinyApp(ui = ui, server = server)
