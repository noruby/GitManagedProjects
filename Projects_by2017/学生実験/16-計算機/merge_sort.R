N <-100000
array <- 1:N
inf <- 1000
data <- cos(-array)
sorted_data <- array(0,dim=c(N))
sorted_data1 <- sort(data)

merge_sort <- function(dat){
	length_dat <- length(dat)
	if(length_dat>1){
		n1 <- round(length_dat/2)
        n2 <- length_dat -n1
        dat1 <- dat[1:n1]
        dat2 <- dat[(n1+1):length_dat]
        dat1 <- merge_sort(dat1)
        dat2 <- merge_sort(dat2)
        dat <- merge(dat1,dat2)
    }
    dat
}

merge <- function(dat1, dat2){
	length1 <- length(dat1)
	length2 <- length(dat2)
	length <- length1+length2
	dat1 <- append(dat1, inf)
	dat2 <- append(dat2, inf)
	da <- rep(0,length)
	i1 <- 1 
	i2 <- 1
    for (i in 1:length){
		if(dat1[i1]<dat2[i2]){
			da[i] <- dat1[i1]
			i1=i1+1
		} else {
			da[i] <- dat2[i2]
			i2=i2+1
		}
	}
	da
}

sorted_data <- merge_sort(data)
xrange <- c(1,N)
yrange <- c(min(data),max(data))
png("merge_sort.png", width = 700, height = 500) 
#plot(data, xlim=xrange, ylim=yrange, col="black")  
#par(new=T)
plot(sorted_data, xlim=xrange, ylim=yrange, col="blue", type="l")  
par(new=T)
plot(sorted_data1, xlim=xrange, ylim=yrange, col="red", type="l")  
dev.off() 
