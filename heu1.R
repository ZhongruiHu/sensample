#----------------------------
# Simulation of heuristics 1
#----------------------------
require('stats')
graphics.off()

#- - - - - - - - BEGIN CONFIGURABLE SECTION - - - - - - - -
datafname <- 'AllDataFromOlga_NellyBay.txt'
sensacc <- 0.3

# alternative data sources:
#datafname <- 'IntelDataSet.txt'
#sensacc <- 0.01

all_data <- scan(datafname)
fin_data <- all_data

beg_idx <- 1
end_idx <- 50
max_idx <- 270

# initialize plot environment
par(mfrow=c(2,4))

# prediction horizon
horiz <- 5

#- - - - - - - - END CONFIGURABLE SECTION - - - - - - - -

tot_skips <- 0
tot_smpls <- 0

z = ts(all_data[beg_idx:end_idx], start=beg_idx)

while(TRUE) {
	plot(z)
	
	#- - - - - - - - Power transformation - - - - - - - -
	# TODO: determine whether to perform power transformation

	#- - - - - - - - Differencing - - - - - - - -
	# TODO: determine d, how many times of differencing to do, for now assume...
	d <- 1
	if(d > 0)
		z = diff(z, differences=d)

	# visual effect of differencing
	plot(z)

	#- - - - - - - - Model identification - - - - - - - -
	acf(z, lag.max=50)
	pacf(z, lag.max=50)

	# TODO: determine p and q, for now, assume...
	# observation: when p=5, the fit is worse, but the prediction is better compared to p=6
	p <- 4
	q <- 0

	#- - - - - - - - Model estimation - - - - - - - -
	model = arima(z, order=c(p,d,q), include.mean=TRUE)
	print(sprintf("sigma=%f", sqrt(model$sigma2)))

	#- - - - - - - - Prediction - - - - - - - -
	fcast <- predict(model, n.ahead=horiz)
	z <- c(z, fcast$pred)						# concatenate predictions to original time series
	xi <- all_data[beg_idx]						# initial value
	z <- diffinv(z, differences=d, xi=xi)		# recover original time series
	z_n_l <- z[(length(z)-horiz+1):length(z)]	# forecasts
	half_intrvl <- 1.96*fcast$se				# half of confidence interval (search 'AirPassengers' in R manual)	
	lb <- z_n_l - half_intrvl					# lower bounds
	ub <- z_n_l + half_intrvl					# upper bounds

	# diagnostics
	for(i in 1:horiz) {
		if(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i]) {
			print(sprintf("(%2d) actual=%f, fcast=%f, err=%f,  in [%f, %f] of size %f",
				i, all_data[end_idx+i], z_n_l[i], abs(all_data[end_idx+i]-z_n_l[i]), lb[i], ub[i], ub[i]-lb[i]))
		} else {
			print(sprintf("(%2d) actual=%f, fcast=%f, err=%f, !in [%f, %f] of size %f",
				i, all_data[end_idx+i], z_n_l[i], abs(all_data[end_idx+i]-z_n_l[i]), lb[i], ub[i], ub[i]-lb[i]))
		}
	}

	# determine how many samples to skip
	skip <- 0
	for (i in 1:horiz) {
		if(half_intrvl[i] <= sensacc)
			skip <- skip + 1
		else
			break
	}
	if(skip > 0)
		print(sprintf("skips %d samples: %d to %d", skip, end_idx+1, end_idx+skip))
	tot_skips <- tot_skips +skip
	tot_smpls <- end_idx

	#- - - - - - - - Resume sampling - - - - - - - -
	# slide the sampling window one sample forward
	if(skip > 0)
		slide_step <- skip
	else
		slide_step <- horiz
	beg_idx <- beg_idx + skip + slide_step
	if(beg_idx > max_idx) break
	z <- all_data[beg_idx:end_idx]
	if(skip > 0) {
		z <- c(z, z_n_l[1:skip])
		fin_data[(end_idx+1):(end_idx+skip)] <- z_n_l[1:skip]
	}
	end_idx <- end_idx + skip + slide_step
	if(end_idx > max_idx) break
	z <- c(z, all_data[(end_idx-1):end_idx])
	z <- ts(z, start=beg_idx)
}

print(sprintf("tot_skips=%d, tot_smpls=%d", tot_skips, tot_smpls))

# compare raw and skip-sampled data
dev.new()
plot(ts(all_data[1:max_idx], start=1), ylab="")
par(col='blue')
lines(ts(fin_data[1:max_idx], start=1))
par(col='black')
