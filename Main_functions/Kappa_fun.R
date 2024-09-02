Kappa_fun_1 <- function(X, bw = 0.001,
                      ox.seq = seq(-1,1, l = 50), 
                      oy.seq = seq(-1,1, l = 50)){
  
  W <- X$window
  eval.grid <- expand.grid(ox.seq, oy.seq)
  eval.grid.ppp <- ppp(x = eval.grid[,1], y = eval.grid[,2], window = W)
  
  r0 <- sqrt( (X$window$xrange[1] - X$window$xrange[2])^2 + (X$window$yrange[1] - X$window$yrange[2])^2)
  
  kappa.loc <- function(i){
    vvv <- mvtnorm::dmvnorm(x = cbind(eval.grid.ppp$x, 
                                       eval.grid.ppp$y), 
                             mean = rbind(X$x[i],X$y[i]), 
                             sigma = diag(bw,2) )
    #vvv <- vvv/mvtnorm::dmvnorm(x = cbind(X$x[i],X$y[i]), 
    #                            mean = rbind(X$x[i],X$y[i]), 
    #                            sigma = diag(bw,2) )
    return(vvv)
  }
  
  
  
  kappa.mat <- rowSums( sapply(seq(1,X$n), kappa.loc) )
  result <- as.im(matrix(kappa.mat, ncol = length(ox.seq), nrow = length(oy.seq), byrow = T),
                  xrange = W$xrange, 
                  yrange = W$yrange)
  return(result)
}



Kappa_fun_2 <- function(X, bw = 0.001,
                      ox.seq = seq(-1,1, l = 50), 
                      oy.seq = seq(-1,1, l = 50)){
  
  W <- X$window
  eval.grid <- expand.grid(ox.seq, oy.seq)
  eval.grid.ppp <- ppp(x = eval.grid[,1], y = eval.grid[,2], window = W)
  
  r0 <- sqrt( (X$window$xrange[1] - X$window$xrange[2])^2 + (X$window$yrange[1] - X$window$yrange[2])^2)
  
  my.cross.dist <- crossdist(eval.grid.ppp, X)
  value <- apply(my.cross.dist, 1, min)
  result <- as.im(matrix(exp(-(value/bw) ), ncol = length(ox.seq), nrow = length(oy.seq), byrow = T),
                  xrange = W$xrange, 
                  yrange = W$yrange)
  return(result)
}

Kappa_fun_3 <- function(X, bw = 0.001,
                      ox.seq = seq(-1,1, l = 50), 
                      oy.seq = seq(-1,1, l = 50)){
  
  W <- X$window
  eval.grid <- expand.grid(ox.seq, oy.seq)
  eval.grid.ppp <- ppp(x = eval.grid[,1], y = eval.grid[,2], window = W)
  
  r0 <- sqrt( (X$window$xrange[1] - X$window$xrange[2])^2 + (X$window$yrange[1] - X$window$yrange[2])^2)
  
  my.cross.dist <- crossdist(eval.grid.ppp, X)
  value <- (my.cross.dist <= bw)
  value <- (apply(value, 1, sum) > 0)*1
  result <- as.im(matrix(value, ncol = length(ox.seq), nrow = length(oy.seq), byrow = T),
                  xrange = W$xrange, 
                  yrange = W$yrange)
  return(result)
}


Kappa_fun <- function(X, bw = 0.001,
                      ox.seq = seq(-1,1, l = 50), 
                      oy.seq = seq(-1,1, l = 50),
                      Kappa_select = Kappa_select){
  V0 <- switch(Kappa_select,
         "Kappa_1" = Kappa_fun_1(X, bw = bw, ox.seq = ox.seq, oy.seq = oy.seq),
         "Kappa_2" = Kappa_fun_2(X, bw = bw, ox.seq = ox.seq, oy.seq = oy.seq),
         "Kappa_3" = Kappa_fun_3(X, bw = bw, ox.seq = ox.seq, oy.seq = oy.seq)
        )
  
  return(V0)
}


# For deaths
Kappa_fun_marks <- function(X, psi = 7,
                            ox.seq = seq(-1,1, l = 50), 
                            oy.seq = seq(-1,1, l = 50)){
  W <- X$window
  eval.grid <- expand.grid(ox.seq, oy.seq)
  eval.grid.ppp <- ppp(x = eval.grid[,1], y = eval.grid[,2], window = W);  #X
  
  #closepairs
  
  my.cross.dist <- crosspairs(eval.grid.ppp, X, rmax = (3*psi), twice = TRUE, what = "ijd")
  my.cross.dist$exp.val <- exp(-(my.cross.dist$d/psi)^2)*X$marks[my.cross.dist$j,2]
  vals <- aggregate(exp.val ~ i, data = my.cross.dist, FUN = sum)
  
  result <- rep(0, length(ox.seq)*length(oy.seq))
  result[vals$i] <- vals$exp.val
   
  result <- as.im(matrix(result, nrow = length(oy.seq), ncol = length(ox.seq), byrow = TRUE),
                  xrange = W$xrange, 
                  yrange = W$yrange)
  return(result)
}

Kappa_fun_marks_min <- function(X, psi = 7,
                            ox.seq = seq(-1,1, l = 50), 
                            oy.seq = seq(-1,1, l = 50)){
  W <- X$window
  eval.grid <- as.matrix( expand.grid(ox.seq, oy.seq) )
  #eval.grid.ppp <- ppp(x = eval.grid[,1], y = eval.grid[,2], window = W)
  
  my.cross.dist <- FNN::get.knnx( data = cbind(X$x, X$y), query = eval.grid, k = 1 )
  my.cross.dist$exp.val <- exp(-(my.cross.dist$nn.dist/(X$marks[my.cross.dist$nn.index,2]*psi))^2)
  #vals <- aggregate(exp.val ~ i, data = my.cross.dist, FUN = sum)
  
  result <- rep(0, length(ox.seq)*length(oy.seq))
  result <- my.cross.dist$exp.val
  
  result <- as.im(matrix(result, nrow = length(oy.seq), ncol = length(ox.seq), byrow = TRUE),
                  xrange = W$xrange, 
                  yrange = W$yrange)
  return(result)
}

# For recruits
Kappa_fun_marks_2 <- function(X, Y, psi = 7,
                            ox.seq = seq(-1,1, l = 50), 
                            oy.seq = seq(-1,1, l = 50)){
  W <- X$window
  
  eval.grid <- expand.grid(ox.seq, oy.seq)
  eval.grid.ppp <- ppp(x = eval.grid[,1], y = eval.grid[,2], window = W)
  
  my.cross.dist <- crosspairs(eval.grid.ppp, X, rmax = (3*psi), twice = TRUE, what = "ijd")
  if( all(!is.null(X$marks)) ){
    my.cross.dist$exp.val <- exp(-(my.cross.dist$d/psi)^2)*X$marks[my.cross.dist$j,2]
  }else{
    my.cross.dist$exp.val <- exp(-(my.cross.dist$d/psi)^2)
  }
  vals <- aggregate(exp.val ~ i, data = my.cross.dist, FUN = sum)
  
  result <- rep(0, length(ox.seq)*length(oy.seq))
  result[vals$i] <- vals$exp.val
  
  result <- as.im(matrix(result, nrow = length(oy.seq), ncol = length(ox.seq), byrow = TRUE),
                  xrange = W$xrange, 
                  yrange = W$yrange)
  return(result)
}

Kappa_fun_marks_sim <- function(X, psi = 7,
                            ox.seq = seq(-1,1, l = 50), 
                            oy.seq = seq(-1,1, l = 50)){
  W <- X$window
  eval.grid <- expand.grid(ox.seq, oy.seq)
  eval.grid.ppp <- ppp(x = eval.grid[,1], y = eval.grid[,2], window = W);  #X
  
  my.cross.dist <- crosspairs(eval.grid.ppp, X, rmax = (3*psi), twice = TRUE, what = "ijd")
  my.cross.dist$exp.val <- exp(-(my.cross.dist$d/psi)^2)*X$marks[my.cross.dist$j]
  vals <- aggregate(exp.val ~ i, data = my.cross.dist, FUN = sum)
  
  result <- rep(0, length(ox.seq)*length(oy.seq))
  result[vals$i] <- vals$exp.val
  
  result <- as.im(matrix(result, nrow = length(oy.seq), ncol = length(ox.seq), byrow = TRUE),
                  xrange = W$xrange, 
                  yrange = W$yrange)
  return(result)
}
