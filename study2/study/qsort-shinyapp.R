  library(shiny) #For Shiny Apps
  library(dragulaR) #For draggable widgets in UI
  library(tidyverse) #For piping code and data wrangling
  library(shinyjs) #For javascript use in Shiny
  library(shinyalert) #For pop-up Shiny messages
  library(RMySQL) #Saving data from app to database
  #ShinyBS for tooltips? look into

  #####
  #Description of script
  #####
  #Script is a Shiny App to execute a Q sort task. Participants read statements and indicate their categorical preference
  #After indicating their preference, participants sort cards into a Q sort grid by dragging and dropping statements
  #Upon completion of the task, participants are asked for their reasoning behind the placement of statements into the extreme categories
  #Script produces a csv detailing the user's categorical preference and placement of each statement, and reasoning provided
  #######


  ###################
  #Inputs provided by user
  ###################

  #Import statements for Q sort
  #Sort by order & extract only the text
   statements <- read_csv("q-statements-ordered.csv", col_names = T,
                    cols(order = col_integer(),
                         id = col_integer(),
                         statement = col_character())) %>%
                 dplyr::arrange(order) %>%
                 pull(statement)
   col.distribution <- c(1,2,4,5,6,5,4,2,1) #distribution of statements in final q.sort, from left to right

  #test sample
  #statements <- c("this", "is", "a", "test", "!")
  #col.distribution <- c(1,1,1,1,1)

  options(mysql = list(
    "host" = "localhost",
    "port" = 3306,
    "user" = "root",
    "password" = "JNkoq*3P67"
  ))
  databaseName <- "qsort"
  table <- "responses"

  #Labels
  label.rating <- "auto" #can choose from "auto" (automatically calculated), "none" (no labels), or enter own (e.g., c("-2", "-1", "0", "+1", "+2"))

  #Graphical
  statement.size <- 70/(max(col.distribution) + 1) #height of each statement (vh)
  statement.size.border <- 0.2 #size of border of each statement widget (vh)

  #Warnings for inconsentiencies
  incon.parameter <- list(threshold = 1, width = 1)

  #Include phase for statement suggestions (for pilot testing)
  istest <- T

  #To change
  #Changes to the scale ("unlike my point of view" to "like my point of view") cannot be altered purely through these inputs
  #Currently, if length(col.distributions) < 4, there'll be issues with the placement of grid headings


  #Participants must explain placement of statements
  #X = horizontal location on grid, i.e. which column. 1 (left) to max (right)
  #Y = vertical location on grid, i.e., where in column. 1 (heighest) to max (lowest)
  statements.explain <-
    data.frame(
      X = c(1, length(col.distribution)),
      Y = c(1, 1)
    )

  #Instructions
  #instructions <- c(
  #  "Welcome to the Q sort ShinyApp! :). This app is a work-in-progress. If you'd like to try out the app, please press the button below to continue.",
  #  "Here is a statement about climate change. Please click the button that best matches your point of view",
  #  "Here is a statement about climate change. Please drag it into the grid below to an area that best matches your point of view. Once you've placed all the statements, a button will appear that allows you to continue",
  #  "Please explain the mess you've just made"
  #)

  instructions <- c(
        'In this task, you will be asked to read and give your opinion on different statements about climate change.<br /><br />Please <strong>maximise your browser window</strong> and <strong>click on the continue button </strong>on the top-left&nbsp;to start the task.<br /><br />You may revisit the instructions at any time by clicking the help button found in the top-left corner.',
        'You will now be shown a series of statements about climate change. For each statement, please click one of three buttons below that best corresponds to your point of view.<br /><br /> Once you have assessed all statements, the continue button will appear. Click the continue button to move to the next stage of the task.',
        'You will now be shown the same statements about climate change. <strong>Drag each statement into the grey grid below, and place it in an area that best corresponds to your point of view. The aim is to order statements from least like your point of view to most like your point of view. If you find you agree with many statements, you may need to place statements you agree with on the left part of your grid.</strong><br /><br />The first statements you will read are ones like your point of view, then you will read statements unlike your point of view, finishing with statements that you are unsure or neutral towards.<br /><br />You may rearrange statements in the grid. If you need space to hold spare cards, place them in the green column on the very right.<br /><br />You will be given the option to continue when all statements are placed onto the gray grid.',
        'Please explain how each of these statements are most like/unlike your point of view.'
        )


  ##Quick checks on inputs

  #Check if col distribution is symmetrical, if not, issue warning
  length(col.distribution) %/% 2 %>%
  {if(any(col.distribution[seq(1, .) + . + 1] != col.distribution[seq(-1, -1 * .) + . + 1])){
    warning("Distribution of Q sort columns are asymmetrical")}
  }

  #Stop code if less than two statements
  if(length(statements) < 2){
    stop("Script requires at least two statements")
  }


  #Stop code if mismatch between number of statements and Q sort grid space
  if(length(statements) != sum(col.distribution)){
    sprintf("Error: The Q sort grid has space for %d statements, whereas you have declared %d statements.
            The spaces in the Q sort grid must equal the number of statements", sum(col.distribution), length(statements)) %>%
      stop()
  }

  #transform statements.explain to suitable format


  ####
  #Functions
  ####

  ####
  #UI-side functions
  ####

  makeElement <- function(content, data, content.colname, id.colname)
  {
    #Creates widgets of each statement, with properties extracted from the statement data frame
    #Content is extracted from the content.colname column
    #Id is extracted from the id.colname column
    #Size determines the height of the div
    div(style = sprintf("border-width:%fvh;border-style:solid; height:%fvh; overflow: hidden;", statement.size.border, statement.size),
        drag = data %>% filter(!!as.symbol(content.colname) == content) %>% select(id.colname) %>% pull,
        id = data %>% filter(!!as.symbol(content.colname) == content) %>% select(id.colname) %>% pull,
        title = content,
        div(class = "active title", content)
    )
  }

  create_col <- function(nameofcol, maxstatements, pixel.per.statement, heading){
    #Creates column for Q sort grid, nameofcol is ID
    #maxstatements is the maximum number of statements possibly allocated to col
    #pixel.per.statement is the height of each slot in each column
    #heading is for individual column headings [currently redudant]
    if(nameofcol != "Free"){
      bg <- "#E0E0E0" #background colour
      w <- 1 #width
    }
    else {
      bg <- "#b4ecb4"
      w <- spare.width.q
    }
    if(nameofcol == "C1"){
      os <- offset.q
    }
    else{
      os <- 0
    }
    column(width = w, offset = os, align = "center",
           #h4(heading),
           div(id = nameofcol,
               style = sprintf("height:%fvh;background-color:%s;padding:0px;margin:0px;",
                               pixel.per.statement*maxstatements, bg)
           )
    )
  }

  ####
  #Server-side functions
  ####

  state_join <- function(state, data){
    #Takes the state of draggable divs and updates the statement data frame
    state %>%
      lapply(as.character) %>%
      stack %>%
      rename(id = values, current = ind) %>%
      left_join(select(data, -current), ., by = 'id') %>%
      return()
  }

  update_statepref <- function(data, level){
    #Update function for categorical preferences (Like, Unlike, Neutral)
    #Input the data df and corresponding level (Like = 1, Unlike = 2, Neutral = 3)
    #Returns updated data df with assigned category
    id_to.update <- data %>%
                    filter(pref == "Nonassigned") %>%
                    slice(1) %>%
                    pull(id)
    data %>%
      mutate(pref = replace(pref, id == id_to.update, levels(data$pref)[level])) %>%
      return()
  }


  calc_structure <- function(statement.data, structure.data){
    #calculates column df based on the current state
    #Firstly calculate current contents of each column
    tally.current <-  statement.data %>%
      select(current) %>%
      count(current) %>%
      mutate(col = as.character(current)) %>%
      select(-current)
    #Join this to the column data, update 'current' col
    structure.data %>%
      mutate(prev = current) %>%
      select(-current) %>%
      left_join(.,tally.current, by = "col") %>%
      rename(current = n) %>%
      return()
  }



  detect.incon <- function(df, w){
    #Takes df.s-like df, and width of detection of inconsistencies, and returns # of incon
    L <- length(col.distribution)
    df %>%
    mutate(incon = (pref == "Like" & current %in% paste("C", seq(1, (w - 1 + 1)), sep = "")) |
             (pref == "Unlike" & current %in% paste("C", seq((L - (w - 1)), L), sep = ""))
    ) %>%
      pull(incon) %>%
      return(.)
  }

  #SQLite function for saving data to database
  # saveData <- function(data) {
  #   # Connect to the database
  #   db <- dbConnect(SQLite(), sqlitePath)
  #   # Construct the update query by looping over the data fields
  #   query <- sprintf(
  #     "INSERT INTO %s (%s) VALUES ('%s')",
  #     table,
  #     paste(names(data), collapse = ", "),
  #     paste(data, collapse = "', '")
  #   )
  #   # Submit the update query and disconnect
  #   dbGetQuery(db, query)
  #   dbDisconnect(db)
  # }

  saveData <- function(data) {
    # Connect to the database
    db <- dbConnect(MySQL(), dbname = databaseName, host = options()$mysql$host,
                    port = options()$mysql$port, user = options()$mysql$user,
                    password = options()$mysql$password)
    # Construct the update query by looping over the data fields
    query <- sprintf(
      "INSERT INTO %s (%s) VALUES ('%s')",
      table,
      paste(names(data), collapse = ", "),
      paste(data, collapse = "', '")
    )
    # Submit the update query and disconnect
    dbGetQuery(db, query)
    dbDisconnect(db)
  }



  #qsconvert function to convert query string ID of 7 digits into another 7 digit sequence
  #fed back to participants to enter back into qualtrics
  source("qsconvert.R")


  ####
  #Initialise data structures for server
  ####

  #Initalise data for statements
  dat.s <-  tibble(
      id = statements %>% length %>% seq(1, .) %>% as.character,
      txt = statements,
      pref = "Nonassigned",
      prev = "Available",
      current = "Available"
    )
  #Declare appropriate variables as factors
  dat.s$pref <- factor(dat.s$pref, levels = c("Like", "Neutral", "Unlike", "Nonassigned"))
  dat.s$prev <- factor(dat.s$prev, levels = c("Available", sapply(1:length(col.distribution), function(x){sprintf("C%d", x)})))
  dat.s$current <- factor(dat.s$current, levels = c("Available", sapply(1:length(col.distribution), function(x){sprintf("C%d", x)})))

  #Initalise data for columns of Qsort grid

  dat.structure <-  tibble(
      col = c("Available", sapply(1:length(col.distribution), function(x){sprintf("C%d", x)}), "Free"), #Available for unallocated statements, 'C' for Q sorted column
      max = c(length(statements), col.distribution, max(col.distribution)), #maximum number of cards for each row
      prev = c(length(statements), rep(0, length(col.distribution) + 1)), #previous allocation in each row
      current = 0 #current allocation in each row
    )

  ####
  #Global variables
  ###

  offset.q <- (12- length(col.distribution) - 1) %/% 2 #left hand offset of qsort grid
  spare.width.q <- 12 - offset.q - length(col.distribution) #length of spare column on right hand side

  ####
  #Initalise UI elements
  ###

  #Create ratings labels
  if(label.rating[1] == "auto"){
    #calculate rating labels, from + X to - X, omit 0 if even number of columns
    l <- length(col.distribution) %/% 2 %>%
          seq(-.,.) %>%
            sapply(function(x){
              if(x > 0){
                return(sprintf("+%d", x))
              }
                return(x)
              })
    if(length(col.distribution) %% 2 == 0){
      l <- subset(l, l != as.character(0))
    }
    print("Calculated rating labels automatically")
    label.rating <- l
  }
  if(label.rating[1] == "empty"){
    print("Omitted rating labels")
    label.rating <- rep("", length(col.distribution))
  }


  #Statement divs
  statement.divs <- lapply(statements, makeElement, data = dat.s, content.colname = 'txt', id.colname = 'id')

  #Column layout
  qsortgrid.div <- dat.structure %>%
    filter(col != "Available") %>%
    pull(col) %>%
    lapply(function(x){
      dat.structure %>%
        filter(col == x) %>%
        pull(max) %>%
        {create_col(x, ., statement.size)}
    })

  #Qsort buttons for preference stage of task
  qsortcategories.div <- dat.s$pref %>%
    levels %>%
    subset(. != "Nonassigned") %>%
    rev %>%
    lapply(function(x){
      buttontext.extreme <- "I would say this statement is"
      buttontext.mid <- "I am <b>neutral or unsure</b> <br> of this statement"
      column(3, align = "center",
             actionButton(sprintf("qsortbutton.%s", x),
                          HTML(if_else(x == "Like",
                                  sprintf("%s <br> <b>like</b> my point of view", buttontext.extreme),
                                  if_else(x == "Unlike",
                                          sprintf("%s <br> <b>unlike</b> my point of view", buttontext.extreme),
                                          buttontext.mid
                                          )
                                  )),
                          style = sprintf("border-width:%fvh;border-style:solid;colour:LightPink;font-size:120%%;align:center", statement.size.border)
                          )
      )
    })


  # insert_proceedtext <- function(stage){
  #   #Text which indicates when users can use the proceed button
  #   msg <- if_else(stage == 3, "Once you are happy with your arragenment of cards, c", "C") %>%
  #     {sprintf("%slick the button in the <b>top right</b> of your browser to proceed.", .)}
  #     insertUI(
  #       selector = "#pane",
  #       where = "beforeBegin",
  #       ui = div(id = "continue", h3(HTML(msg)))
  #       )
  # }
  #

  #####
  #Create app
  ####


  ui <- fluidPage(
    div(id = "all",
    shinyjs::useShinyjs(),
    useDragulajs(),
    useShinyalert(),  # Set up shinyalert

    tags$head(tags$style(HTML("#all {
                             display: none;
    }"))),

    #, #qsort.ratingscale, #qsort.ratingscale-num, #qsort.grid, #qsort.buttons

    fluidRow(style = "margin: 0px;", id = "qsort.instructions",
             column(6, align = "left",
                    actionButton("help", h5("Help"))
             ), #Button to re-read instructions
             column(6, align = "right",
                    actionButton("continue", h5(htmlOutput("continue.text")))
                    ) #Button for completing Q sort stages
    ),
    #column(11,
    fluidRow(style = "margin: 0px;padding:0px",
             column(4, offset = 4, align = "center",
                    div(id = "pane")),
                    # div(id = "s",
                    #     style = sprintf("height:%fvh;padding:0px;border-style: inset;", (statement.size+2*statement.size.border)),
                    #     statement.divs[[1]]),
                    # fluidRow(style = "margin: 0px;padding:0px",
                    #          column(9,
             column(3, offset = 4, align = "center",
                    div(id = "Available",
                        style = sprintf("height:%fvh;padding:0px;border-style: inset;", (statement.size+2*statement.size.border)),
                        statement.divs)
             )#) #) #creates a reading pane to view statements
    ),

    fluidRow(style = "margin: 0px;padding:0px; align:center", id = "qsort.ratingscale",
             column(2, offset = offset.q,
                    align = "left",
                    h4(HTML("Statements <b>most unlike</b> my point of view"))),
             column(2, offset = (length(col.distribution) - (2 * 2)), align = "right",
                    h4(HTML("Statements <b>most like</b> my point of view"))),
             column(width = spare.width.q,
                    offset = 0,
                    align = "center",
                    h4("Spare statements")
                    )
             ), #headings for Q sort grid columns
    fluidRow(style = "margin: 0px;padding:0px; align:center", id = "qsort.ratingscale-num",
             column(offset.q),
             label.rating %>%
               lapply(., function(x){
                 column(1,
                        align = "center",
                        h5(as.character(x))
                        )
               })

    ), #rating scale numbers for Q sort grid

    fluidRow(style = "margin: 0px;padding:0px; align:center", id = "qsort.grid",
      qsortgrid.div
    ), #the grid itself, to which statements can be dragged

    fluidRow(style = "margin:0px;padding:100px", id = "qsort.buttons",
            column(1),
            qsortcategories.div
    )#, #categorical sorting buttons for indicating preference

    #Debug central
    # fluidRow(column(12,
    #                 tableOutput('table2')
    #          )
    # ),
    #
    # #Debug central
    # fluidRow(column(12,
    #                 tableOutput('table1')
    # )
    # ),
    ,

    dragulaOutput("dragula"),
    extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }", functions = c("closeWindow"))
  ))

  server <- function(input, output, session) {
    #Declared variables
    qsort.stage <- reactiveVal(1) #Stage of task
    df.s <- reactiveVal(dat.s) #data for each statement
    df.structure <- reactiveVal(dat.structure) #data for each column in the q sort grid
    rv <- reactiveValues(expl = list(most.like = NULL, most.unlike = NULL)) #open-ended questions following Q sort
    wv <- reactiveValues(noroom = 0, incon = 0, incon.error = 0) #list of warning counts
    supress <- reactiveValues(noroom = F, incon = F, finalsort = F) #list of supression options, for room in grid and inconsistency
    qsort.time <- reactiveValues(time = rep(Sys.time(),6)) #holds times each qsort stage began
    pilot.ans <- reactiveValues(p1 = "", p2 = "") #responses for pilot questions
    place <- reactiveValues(avail = 0, grid = 0,  spare = 0, spare.max = 0) #for monitoring use of area where statements are placed, keep track of max. statements in spare statement area at any one time
    current <- reactiveValues(noroom = F, begunsort = F) #current states of UI
    qs <- reactiveValues(full = list(0), id = "", pass = "", error = T) #query string

    #Debug
    # output$table1 <- renderTable(df.s())
    # output$table2 <- renderTable(
    #   # tibble(
    #   #   grid = place$grid,
    #   #   avail = place$avail,
    #   #   spare = place$spare,
    #   #   max = place$spare.max
    #   # )
    #   tibble(
    #     id = qs$id,
    #     pass = qs$pass,
    #     error = qs$error
    #   )
    # )

    #Creates a dragable interface
    output$dragula <- renderDragula({
      dragula(dat.structure$col)
    })

    #Reactive variable which is easier to read than input$dragula
    state <- reactive({
      dragulaValue(input$dragula)
    })

    #Hide most of the UI, bar title, intrusctions, and continue button
    phase1UI <- c("help", "continue") #UI for phase 1, not hidden on phase 2
    phase2UI <- c("qsortbutton.Like", "qsortbutton.Unlike", "qsortbutton.Neutral", "Available") #UI for phase 2
    phase3UI <- c("qsort.ratingscale", "qsort.grid", #"qsort.readingpane", "qsort.ratingscale-num") #UI for phase 3
                  "Available", "qsort.ratingscale-num")
    phase4UI <- c("qsort.explain1", "qsort.explain2")
    phase5UI <- c("pilot1", "pilot2")

    #Check qsort.stage() and update UI accordingly
    observeEvent(qsort.stage(), {
      qsort.time$time[qsort.stage()] <- Sys.time()

      #show instructions
      if(qsort.stage() <= length(instructions)){
      shinyalert(title = "Instructions",
                 text = instructions[[qsort.stage()]],
                 html = T,
                 animation = F
                 )
      }


      if(qsort.stage() == 1){
        #hide elements
        c(phase1UI, phase2UI, phase3UI, "end") %>%
          lapply(shinyjs::hide)
        #Retrieve query string, error if length incorrect
        qs$full <- getQueryString()
        qs$error <- qs$full %>%
          {(length(.) != 1)}
        if(!qs$error){
          qs$id <- qs$full$ID
          qs$error <- qs$id %>%
          {(nchar(.) != 7)}
        }
        if(!qs$error){
          qs$pass <- qsconvert(qs$id)
        }
        #Display relevant elements
        shinyjs::show("all")
        if(qs$error){
          insertUI(
            selector = "#pane",
            where = "afterBegin",
            ui = div(id = "error", h1(HTML("<b><u>Error</u></b>")), h2("You have entered the wrong website. This is likely caused by copying and pasting. Please return to Qualtrics and copy and paste the entire link (and only the link)."))
          )
        }
        else {
          phase1UI %>%
            lapply(shinyjs::show)
        }
      }



      if(qsort.stage() == 2){
        shinyjs::hide(id = "continue")
        phase2UI %>%
          lapply(shinyjs::show)
        #Hide all statements except the first
        isolate(df.s()) %>%
          slice(-1) %>%
          pull(id) %>%
          lapply(shinyjs::hideElement)
      }
      if(qsort.stage() == 3){
        #Remove phase2 buttons to prevent whitespace between reading pane and grid
        phase2UI[phase2UI != "Available"] %>%
        # phase2UI %>%
          lapply(function(x){
            y <- gsub(".", "\\.", x, fixed=TRUE) #Jquery uses '.', need to replace with '\\.
            removeUI(selector = sprintf("#%s", y), immediate = T)}) #remove buttons
        hide(id = "continue")
        # insertUI(selector = "#pane",
        #          where = "afterBegin",
        #          ui = create_statement_divs(3))
        phase3UI %>%
          lapply(shinyjs::show)
      }
      if(qsort.stage() == 4){
        #Create post q.sort questions
        insertUI(
          selector = "#pane",
          where = "afterBegin",
          ui = div(list(textAreaInput("qsort.explain1",
                                  label =
                                    sprintf("Please explain how the statement \"<i>%s</i>\" is most like your point of view",
                                            filter(df.s(), df.structure()$col[nrow(df.structure())-1] == current) %>%
                                              pull(txt)) %>%
                                    HTML(),
                                  width = '100%', height = '100%'),
                    textAreaInput("qsort.explain2",
                                  label =
                                    sprintf("Please explain how the statement \"<i>%s</i>\" is most unlike your point of view",
                                            filter(df.s(), df.structure()$col[2] == current) %>%
                                              pull(txt)) %>%
                                    HTML(),
                                  width = '100%', height = '100%')
          ))
        )
        #Remove phase3 buttons
        phase3UI %>%
          lapply(function(x){
            y <- gsub(".", "\\.", x, fixed=TRUE) #Jquery uses '.', need to replace with '\\.
            removeUI(selector = sprintf("#%s", y), immediate = T)}) #remove buttons
      }
      if(qsort.stage() == 5){
       if(istest){
         insertUI(
           selector = "#pane",
           where = "afterBegin",
           ui = list(textAreaInput("pilot1",
                                   label = "In this task, you read and rated a series of statements. Are there any statements crucial to your perspective on climate change which you think we missed?",
                                   width = '100%', height = '100%'
                                   ),
                     textAreaInput("pilot2",
                                   label = "How did you find this task? Do you have any recommendations for improving the experience?",
                                   width = '100%', height = '100%'
                     )
         )
         )
        #Remove phase4 UI
         phase4UI %>%
           lapply(function(x){
              y <- gsub(".", "\\.", x, fixed=TRUE) #Jquery uses '.', need to replace with '\\.
              removeUI(selector = sprintf("#%s", y), immediate = T) #remove buttons
              removeUI(selector = sprintf("label[for='%s']", y), immediate = T) #remove labels
           })
       }
        else{
          #skip
          qsort.stage(qsort.stage() + 1)
        }
      }
      if(qsort.stage() == 6){
        #Task is complete, save results and close app
        #Hide elements
        c(phase1UI, phase5UI) %>%
          lapply(shinyjs::hide)
        #Show passcode
        insertUI(
          selector = "#pane",
          where = "afterBegin",
          ui = div(id = "complete", h1(HTML(sprintf("Your code is: &quot;<b>%s<b>&quot;", qs$pass))), h2("You have completed the task. Please return to Qualtrics and enter your code into the text box (do not include the quotation marks)."))
        )


        #save statement information
        results <- df.s() %>%
          select(-txt, -prev) %>%
          mutate_if(is.factor, as.character) %>%
          rename(sort = current) %>%
          mutate(reason = ifelse(sort == "C1", rv$expl$most.unlike,
                                 ifelse(sort == dat.s$current %>% levels() %>% tail(1), rv$expl$most.like, NA))) %>%
          gather(phase, value, -id) %>%
          mutate(id = sprintf("%s_sta_%s", phase, id)) %>%
          select(-phase) %>%
          spread(id, value)
        #save timing data
        timing <- tibble(
          phase = c("instruct", "pref", "sort", "reason", "pilot"),
          time = sapply(1:5, function(x){difftime(qsort.time$time[x+1], qsort.time$time[x], units = "secs")})
        ) %>%
          spread(phase, time) %>%
          rename_all(funs({sprintf("time_%s", .)}))
        #save functionality data
        #includes errors, and use of spare statement secion
        funct <- tibble(
          count_noroom = wv$noroom,
          count_incon = wv$incon.error,
          supress_noroom = supress$noroom,
          supress_incon = supress$incon
        )
        #save pilot data, label NA if not applicable
        if(istest){
        pilot <- tibble(
          pilot_response_1 = pilot.ans$p1,
          pilot_response_2 = pilot.ans$p2
        )
        }
        else{
          pilot <- tibble(
            pilot_response_1 = NA,
            pilot_response_2 = NA
          )
        }
        #save statement placement data, for future improvement of app
        place <- tibble(
          place_grid = place$grid,
          place_avail = place$avail,
          place_spare = place$spare,
          place_spare_max = place$spare.max
        )
        #save data required to link Qualtrics and ShinyApp
        qualtrics <- tibble(
          id = qs$id,
          passcode = qs$pass,
          time_start = as.character(qsort.time$time[1])
        )
        #save
        #MySQL saveData function does not work if text data contains an apostrophe, convert to &rsquo; (HTML)
        sub_apostrophe <- function(x){return(gsub('\'', '&rsquo;' , x))}
        results %>%
          bind_cols(timing) %>%
          bind_cols(funct) %>%
          bind_cols(pilot) %>%
          bind_cols(place) %>%
          bind_cols(qualtrics) %>%
          mutate_at(vars(starts_with("reason"), starts_with("pilot_response")), sub_apostrophe) %>%
          #write_delim('results.csv', delim = ',')
          saveData()
        #stop the application
        #stopApp()
      }
    })

    #Instructions
    observeEvent(input$help, {
      #show instructions
      if(qsort.stage() < 6){
        shinyalert(title = "Instructions",
                   text = instructions[[qsort.stage()]],
                   html = T,
                   animation = F
        )
      }
    })
    #Continue button text
    output$continue.text <- renderUI({
      s <- if_else(istest, 5, 4) #final stage for user interaction
      if_else(qsort.stage() < s, "Click here to continue", "Click here to complete the study") %>%
        {sprintf("<b>%s</b>", .)} %>%
        HTML()
    })

    ####
    #Server-side functions for qsort stage 1 [Introductory instructions]
    ###

    #No functions as of yet

    ###
    #Server-side functions for qsort stage 2
    ###

    #The following three observers watch the button input for preference and update df.s() accordingly
    observeEvent(input$qsortbutton.Like, {
      df.s() %>%
        update_statepref(1) %>%
        df.s(.)
    }, ignoreInit = T)

    observeEvent(input$qsortbutton.Unlike, {
      df.s() %>%
        update_statepref(3) %>%
        df.s(.)
    }, ignoreInit = T)

    observeEvent(input$qsortbutton.Neutral, {
      df.s() %>%
        update_statepref(2) %>%
        df.s(.)
    }, ignoreInit = T)

    observeEvent(
      #Update display of available statements to the next nonassigned statement
      c(
        input$qsortbutton.Like,
        input$qsortbutton.Unlike,
        input$qsortbutton.Neutral
      ), {
        #Hide all statements
        df.s() %>%
          pull(id) %>%
          lapply(shinyjs::hideElement)
        #Show next unassigned statement
        df.s() %>%
          filter(pref == "Nonassigned") %>%
          slice(1) %>%
          pull(id) %>%
          lapply(showElement)
        #If no unassigned statements, allow user to continue and hide redunant UI
        if(!("Nonassigned" %in% df.s()$pref)){
          phase2UI %>%
            lapply(shinyjs::hide)
          shinyjs::show(id = "continue")

        }
      }, ignoreInit = T)

    ####
    #Server-side functions for qsort stage 3 [Drag into grid]
    ####

    observe({
      #Observes state, and updates df.s() accordingly
      if(qsort.stage() == 3){
        if(!(is_empty(unlist(state())))){
          isolate(df.s()) %>%
            mutate(prev = current) %>%
            state_join(state(), .) %>%
            df.s()
        }
        #check if qsort has begun
        if(!isolate(current$begunsort)){
          current$begunsort <- isolate(df.s()) %>%
            filter(current != "Available") %>%
            nrow() > 0
        }
      }
    })

    observe({
      #When df.s() is updated, updates the reading pane
      if(qsort.stage() == 3){
        next.available <-  df.s() %>%
                              mutate(pref = fct_relevel(pref, c('Like', 'Unlike', 'Neutral'))) %>%
                                filter(current == 'Available') %>% #Reorder the pane, present Liked items first, Unliked second, Neutral last
                                arrange(pref, current) %>%
                                slice(1) %>%
                                pull(id)
        showElement(next.available) #Show the next statement in reading pane
        df.s() %>%
          filter(current == 'Available') %>%
          filter(id != next.available) %>%
          pull(id) %>%
          lapply(shinyjs::hideElement) #Hide all other unallocated statements
        #If all statements are placed, show the end button
        if(!("Available" %in% df.s()$current) & !("Free" %in% df.s()$current)){
          shinyjs::show(id = "continue")
        }
        else{
          shinyjs::hide(id = "continue")
        }
      #Check for inconsistencies between stage 2 and 3 of qsort
      wv$incon <- df.s() %>%
        detect.incon(incon.parameter$width) %>%
        {which(.)} %>%
        length #number of inconstencies
      }
    })

    observe({
      if(qsort.stage() == 3){
        #Update column data
        isolate(df.structure()) %>%
          calc_structure(df.s(), .) %>%
          df.structure()
        #Find any errors, where number of statements in a column exceeeds the maximum
        error <- isolate(df.structure()) %>%
          mutate(overmax = current > max) %>%
          filter(overmax) %>%
          pull(col)

        if(!(is_empty(error))){
          #Then one of the cols > max, need to revert back to previous
          culprit <- isolate(df.s()) %>%
            filter(prev != current)
          #Remove the offending statement from the UI
          removeUI(selector = sprintf("#%s", pull(culprit, id)), immediate = T)
          location <- if_else(
            pull(culprit, prev) == "Available",
            "afterBegin",
            "beforeEnd"
          )
          #Re-insert the offending statement in the previous position
          insertUI(
            selector = sprintf("#%s", pull(culprit, prev)),
            where = location,
            ui = culprit %>%
              pull(id) %>%
              as.integer %>%
              statement.divs[.],
            immediate = T
          )
          #Refresh the grid to update state()
          js$refreshDragulaR("dragula")
          #Error message pops up, if not supressed
          if(!(supress$noroom)){
          if_else(isolate(wv$noroom) < 3, F, T) %>% #show cancel button on error message to surpress further warnings?
          {shinyalert("",
                      paste(c("There is no room in this column to place another statement. You can make room by moving the statements already placed in the grid.", "If you wish to hide this message in the future, press the 'Cancel' button")[c(T, .)], sep = "", collapse = ""),
                      type = "error",
                      closeOnClickOutside = F,
                      closeOnEsc = F,
                      showCancelButton = .,
                      animation = F,
                      callbackR = function(x){supress$noroom <- !(x)}
          )}
            wv$noroom <- isolate(wv$noroom) + 1
          }
          current$noroom <- T #next iteration is to fix noroom error
        }
        else{
          #if noroom error remains untriggered, consider this a valid move
          #ignore if this move is a result of noroom error
          #only begin counting when Qsorting has begun
          if(!isolate(current$noroom) & isolate(current$begunsort)){
            #record user interaction with grid/spare statement area
            #was this statement moved to the spare statement area?
            movedtofree <- isolate(df.s()) %>%
              transmute(prev != "Free" & current == "Free") %>%
              pull() %>%
              any()
            if(movedtofree){
              #incriment count
              place$spare <- isolate(place$spare) + 1
              #note number of statements in spare statement area
              current.free <- isolate(df.s()) %>%
                filter(current == "Free") %>%
                nrow
              #update max if required
              if(isolate(place$spare.max) < current.free){
                place$spare.max <- current.free
              }
            }
            #check if interaction moved to "Available"
            movedtoavail <- isolate(df.s()) %>%
              transmute(prev != "Available" & current == "Available") %>%
              pull() %>%
              any()
            if(movedtoavail){
              #incriment count
              place$avail <- isolate(place$avail) + 1
            }
            if(!movedtoavail & !movedtofree){
              #was moved to grid, incriment count
              place$grid <- isolate(place$grid) + 1
             }
          }
          current$noroom <- F #reset noroom detector
        }
      }
    })

    observeEvent(wv$incon, {
      #If error is relevant
      if(qsort.stage() == 3 & (wv$incon >= incon.parameter$threshold)){
        #AND if inconsistency is new
        new.incon <- which(detect.incon(filter(isolate(df.s()), prev != current), 1) == T) %>%
          length()
        if(new.incon > 0 & !current$noroom){
          wv$incon.error <- isolate(wv$incon.error) + 1
          #Give error msg
          if(!(supress$incon)){
            shinyalert("",
                       "There is an inconsistency between your previous ratings of these statements and where you are placing the statements in the grid. <b>If this is intentional, please ignore this warning</b>.<br><br> You can hide this message in the future by pressing the 'Cancel' button",
                       closeOnClickOutside = F,
                       closeOnEsc = F,
                       showCancelButton = T,
                       animation = F,
                       html = T,
                       callbackR = function(x){supress$incon <- !(x)}
            )
          }
        }
      }
    })

    observe({
      if(!supress$finalsort){
        error <- isolate(df.structure()) %>%
          mutate(overmax = current > max) %>%
          filter(overmax) %>%
          pull(col)
        if(!("Available" %in% df.s()$current) & !("Free" %in% df.s()$current) & ((is_empty(error)))){
          #If first time participant has made complete grid + noroom error has not triggered, give P info
          shinyalert(title = "",
                     text = 'Now you have placed all statements on the grid. Please go over your arrangement once more and shift statements if you would like to. Once you are happy with your arrange, click the continue button',
                     animation = F
                     )
          supress$finalsort <- T #Supress this info in future
        }
      }
    })


    ####
    #Server-side functions for qsort stage 4 [Explanatory questions]
    ###

    #Stage ends when continue button is clicked
    observeEvent(input$continue, {
      if(qsort.stage() < 4){
        advance <- T #Can advance to next stage
      }
      else{
        if(qsort.stage() == 4){
        #Set explanation variables
        rv$expl$most.like <- input$qsort.explain1
        rv$expl$most.unlike <- input$qsort.explain2
        #Answered all questions?
        answered.all <- rv$expl %>%
          unlist %>%
          {any(. == "")} %>%
          {!(.)}
        advance <- answered.all
        if(!answered.all){
          #Can only advance if all questions are answered
          shinyalert("", "Please answer all questions", type = "error", animation = F)
          }
        }
        if(qsort.stage() == 5){
          if(istest){
          #Record pilot answers if responses exist
          pilot.ans$p1 <- input$pilot1
          pilot.ans$p2 <- input$pilot2
          }
          advance <- T #advance to next stage
        }
      }
      if(advance){
        qsort.stage(qsort.stage() + 1)
      }
    })

  }

  shinyApp(ui = ui, server = server, options = c(port = 3838) )
