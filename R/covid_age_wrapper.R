library(data.table)

run_covid_age <- function(covid_model, par_list = NULL, delete_output = TRUE){
  
  ## covid_model: full path to the executable
  ## par_list: a named list that will be passed to the simulation
  
  if(!is.null(par_list)){
    ## parse parameter list
    npar <- length(par_list)
    
    par_names <- names(par_list)
    par_vals <- unlist(par_list)
    
    ## lay out name:val for system command
    input_char <- rep(NA, npar)
    for (ii in 1:npar){
      if(par_names[ii] %in% c("output_directory", "output_filename")){
        input_char[ii] <- paste0('"', par_names[ii], '"', ':"' , par_vals[ii], '"')
      } else {
        input_char[ii] <- paste0('"', par_names[ii], '"', ":" , par_vals[ii])
      }
    }
    input_str <- paste(input_char, collapse = ",")
    input_str <- paste0("'{", input_str , "}'")
    
    ## check if output_file and output_directory are present
    if("output_directory" %in% par_names){
      if("output_filename" %in% par_names){
        out_file <- paste0(par_list["output_directory"], "/", par_list["output_filename"])
      } else {
        out_file <- paste0(par_list["output_directory"], "/", "daily_output.txt")
      }
    } else {
      if("output_filename" %in% par_names){
        out_file <- paste0(getwd(), "/", par_list["output_filename"])
      } else {
        out_file <- paste0(getwd(), "/", "daily_output.txt")
      }
    }
  }
  
  x <- paste(covid_model, input_str)
  # cat(x)
  system(x)
  # cat(out_file)
  ## read the  output
  dt <- fread(out_file)
  
  if(delete_output)  unlink(out_file)
  return(dt)
}

# covid_model <- "/home/afadikar/work/projects/git/covid-age/exp/chicago_yr1/model"
# par_list <- list("random_seed" = 23,
#                  "ini_Ki" = 1.5)
# run_covid_age(covid_model, par_list)
