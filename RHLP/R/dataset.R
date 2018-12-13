MyData <- setRefClass(
  # Set the name for the class
  "MyData",

  # Define the fields
  fields = list(
    x="matrix",
    y="matrix",
    m="numeric"),

  # Set the methods
  methods=list(
    # set Data from a file
    setDataFromMat = function(fileName){
      library(R.matlab)
      data <- readMat(fileName)
      x <<- data$x
      y <<- data$y

      if (ncol(y)!=1){
        y<<-t(y)
      }
      setDataProperties()
    },


    setDataProperties = function(){
      m <<- length(y)
    }
  )
)
