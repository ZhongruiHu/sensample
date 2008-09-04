require('stats')

#- - - - - - - - BEGIN CONFIGURABLE SECTION - - - - - - - -

datafname <- 'AllDataFromOlga_NellyBay.txt'
#datafname <- 'IntelDataSet.txt'

all_data <- scan(datafname)

beg_idx <- 1
end_idx <- 100

par(mfrow=c(3,2))

max_p <- 6
max_q <- 6

#- - - - - - - - END CONFIGURABLE SECTION - - - - - - - -

z <- ts(all_data[beg_idx:end_idx], start=beg_idx)
z0 <- z
plot(z)

d <- 2
if(d > 0)
	z <- diff(z, differences=d)
plot(z)

acf(z, lag.max=50)
pacf(z, lag.max=50)

min_aic <- 1e+6
for(p in 0:max_p) {
	for(q in 0:max_q) {
		if(p == 0 && q == 0) next
		model <- NULL
		tryCatch(model <- arima(z, order=c(p,d,q), include.mean=TRUE),
			error=function(e) {
				writeLines(sprintf("[p=%d, q=%d] ", p, q), sep="");
				writeLines(toString(e), sep="")
			})
		if(!is.null(model)) {			
			write(sprintf("[p=%d, q=%d] aic=%f, sigma=%f", p, q, model$aic, sqrt(model$sigma2)), "")
			if(model$aic < min_aic) {
				min_aic <- model$aic
				best_model <- model
			}
		}
	}
}
writeLines("Best model: ", sep="")
print(best_model)

z <- diffinv(z, differences=d, xi=c(z0[1:d]))
plot(z0, ylab="Reconstructed z")
par(col='blue')
lines(z)
par(col='black')
