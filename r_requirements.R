installed <- installed.packages()[, 1]
required <- c('shiny'
              , 'shinythemes'
              , 'shinyWidgets'
              , 'data.table'
              , 'ggplot2'
              , 'ggthemes'
              , 'gridExtra'
              , 'reticulate')

for (package in required){
  if (!(package%in% installed)){
    install.packages(package)
  }
}
