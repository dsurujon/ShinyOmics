get_ci <- function(x){
  error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
  ci95 <- data.frame('y'=mean(x),'ymin'=mean(x)-error,'ymax'=mean(x)+error)
  return(ci95)
}