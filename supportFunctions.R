## Support functions

# define negation of %in% operator, returns logical vector for 'NOT IN'
'%ni%'=Negate('%in%')

# allow natural sorting on multiple columns
multiMixedOrder=function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}