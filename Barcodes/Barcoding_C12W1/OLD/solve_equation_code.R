
# set-up a test data.frame
test.data <- data.frame(no_observations = c(2, 12, 12, 4),
                        allele = c("A1", "A2", "B1", "B2"),
                        equation = c("A1",
                                     "A2 * (A1 / (A1 + B2))",
                                     "B1 * (B2 / (A1 + B2))",
                                     "B2") )
View(test.data)

# write a function to solve the equation and add a column of solutions to the data.frame

# args

# data - dataset in the format that Bjorn sent i.e. a column for n-observations, allele name and equation
# n_obs - column name containing the number of observations
# allele_col - column name containing the allele names
# equation - column name containing the equations to be solved

solve_allele_equation <- function(data, n_obs, allele_col, equation) {
  
  # create a data.frame tol work on
  df <- data
  
  # assign the number of obsrevations to the correct allele
  for (j in 1:nrow(df)) {
    assign(x = df[[allele_col]][j], value = df[[n_obs]][j] )
  }
  
  # solve the equations
  x <- sapply(df[[equation]], function(x) eval(parse(text = x)) )
  names(x) <- NULL
  
  # add the equation solutions to the original data.frame
  df[["equation_solution"]] <- x
  
  return(df)
  
}

# test the function on the test.data
solve_allele_equation(data = test.data, 
                      n_obs = "no_observations", 
                      allele_col = "allele", 
                      equation = "equation")

### END
