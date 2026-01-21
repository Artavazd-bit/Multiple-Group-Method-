foo <- function(){
  
  #warning("Obbacht!")
  log(-1)
  sqrt(-1)
  #stop("test")
  #warning("Gefahr!")
 # stop("error")
}

foo2 <- function(){
  warning("Obbacht!2")
  warning("Gefahr!2")
  2
}

testfunc()

errormsg <- NA
res <- tryCatch( foo(), 
          warning = function(w) {
            print(conditionMessage(w))
            }, 
          error = function(e){
            print(conditionMessage(e))
            },
          finally = {print("finally reached")}
          ) 

w <- list()
withCallingHandlers(foo(), warning = function(w){w <<- list(w, list(w)); invokeRestart("muffleWarning") })
str(w)
rm(w)



trycatchmore <- function(statfunc){
  warnings <- list()
    tryCatch(
      {
        message("ill try to run the function foo")
        statfunc()
      },
      error = function(cond){
        #message("error occured:")
        #error <- conditionMessage(cond)
        #print(error)
        return(3)
      },
      warning = function(cond) {
        warnings <<- c(warnings, list(cond))
        #message("warning occured:")
        #warning <- conditionMessage(cond)
        return(list(4, warnings) )  
      },
      finally = {
        #message("ich write nonetheless")
      }
  )
}

testcall <- function(statfunc){
  warnings <- list()
  errors <- list()
  withCallingHandlers(
    {
      statfunc()
    },
    error = function(cond){
      errors <<- c(errors, list(cond))
      return(errors)
    },
    warning = function(cond){
      warnings <<- c(warnings, list(cond))
      return(warnings)
    }
  )
}
