require('stats')

#- - - - - - - - BEGIN CONFIGURABLE SECTION - - - - - - - -

if(!exists("which_data"))
	stop('Must specify which_data')
switch(which_data,
	# 1
	{data_fname <- 'Olga'; data_type <- 'Temperature'; err_tol <- 0.3},
	# 2
	{data_fname <- 'Intel1'; data_type <- 'Temperature'; err_tol <- 0.3},
	# 3
	{data_fname <- 'Intel2'; data_type <- 'Humidity'; err_tol <- 0.3},
	# 4
	{data_fname <- 'Intel3'; data_type <- 'Light'; err_tol <- 3},
	# 5
	{data_fname <- '41001h2007WSPD'; data_type <- 'Wind speed';	err_tol <- 0.3},
	# 6
	{data_fname <- '41001h2007GST'; data_type <- 'Gust speed'; err_tol <- 0.3},
	# 7
	{data_fname <- '41001h2007WVHT'; data_type <- 'Wave height'; err_tol <- 0.3},
	# 8
	{data_fname <- '41001h2007PRES'; data_type <- 'Pressure'; err_tol <- 1})
all_data <- scan(data_fname)

beg_idx <- 1
end_idx <- 1000 # 0 = till the end

par(mfrow=c(3,2))

max_p <- 6
max_q <- 6

#- - - - - - - - END CONFIGURABLE SECTION - - - - - - - -

if(end_idx > 0) {
	z <- ts(all_data[beg_idx:end_idx], start=beg_idx)
} else {
	z <- ts(all_data[-(1:(beg_idx-1))], start=beg_idx)
}
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
