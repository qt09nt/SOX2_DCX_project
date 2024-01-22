### In this file I replicate the functions present in perseus based on the description in coxdocs


# implements the method of width adjustment from perseus on the matrix "data"
# returns modified matrix, does not modify the original
# (http://www.coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:normalization:widthadjustment)
width_adjustment = function(data_in){
  data = as.matrix(data_in)
  q = quantile(data, na.rm = TRUE)[2:4] # only care about indices 2:4 because we don't need min/max
  data = data - q[2]
  data[data > 0 && !is.na(data)] = data[data > 0 && !is.na(data)] / (q[3] - q[2])
  data[data < 0 && !is.na(data)] = data[data < 0 && !is.na(data)] / (q[2] - q[1])
  return(data)
}

z_score_by_row = function(data_in){
  data = as.matrix(data_in)
  for(row in 1:nrow(data)){
    m = mean(data[row,], na.rm=T)
    s = sd(data[row,], na.rm=T)
    data[row,] = data[row,] - m
    data[row,] = data[row,] / s
  }
  return(data)
}

z_score_by_col = function(data_in){
  data = as.matrix(data_in)
  for(col in 1:ncol(data)){
    m = mean(data[,col], na.rm=T)
    s = sd(data[,col], na.rm=T)
    data[,col] = data[,col] - m
    data[,col] = data[,col] / s
  }
  return(data)
}

# filter rows that have a minimum percent valid (non-NA values)
filter_valid = function(data, percent = 100){
  num_valid = rowSums(!is.na(data))
  num_cols = dim(data)[2]
  rows_percent_valid = num_valid / num_cols
  data_out = data[rows_percent_valid >= (percent / 100),]
  return(data_out)
}


select_columns = function(data, colname1, colname2){
  first_col = which(colnames(data) == colname1)
  last_col = which(colnames(data) == colname2)
  return(data[,first_col : last_col])
}


# assumes length(vector1) == length(vector2)
scatter_plot = function(vector1 , vector2){
  data = cbind.data.frame(d1 = vector1, d2 = vector2)
  p = ggplot(data = data, mapping = aes(x = d1, y = d2)) + geom_point()
  return(p)
}

# implements the "Replace missing values from normal distribution" function in perseus
# http://www.coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:imputation:replacemissingfromgaussian
repl_NA_from_norm = function(data, width = 0.3, down_shift = 1.8,
                             mode = c("columns", "whole")){
  
  if(mode == "whole"){
    stdev = sd(as.matrix(data), na.rm = TRUE) * width
    mean = mean(as.matrix(data), na.rm = TRUE) - (down_shift * stdev)
    num_values = sum(is.na(data))
    generated_values = rtnorm(num_values, mean, stdev)
    data_out = data
    data_out[is.na(data)] = generated_values
  }else if(mode == "columns"){
    data_out = data
    for(i in 1:ncol(data)){
      colsd = sd(as.matrix(data[,i]), na.rm = TRUE) * width
      colmean = mean(as.matrix(data[,i]), na.rm = TRUE) - (down_shift * colsd)
      num_values = sum(is.na(data[,i]))
      generated_values = rtnorm(num_values, colmean, colsd)
      data_out[is.na(data[,i]),i] = generated_values
    }
  }else{
    stop("mode = c(\"columns\", \"whole\")")
  }
  return(data_out)
}


