# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation of various heuristics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
require('stats')
graphics.off()

#- - - - - - - - BEGIN CONFIGURABLE SECTION - - - - - - - -

datafname <- 'Olga'
datatype <- 'Temperature'
errtol <- 0.3
#datafname <- 'Intel'
#datatype <- 'Temperature'
#errtol <- 0.5

all_data <- scan(datafname)
fin_data <- all_data

beg_idx <- 1
end_idx <- 50
max_idx <- 100

# maximum p and q
max_p <- 4
max_q <- 4

# forecast horizon
horiz <- 5

# forecast limit
fcast_lim <- 0.90
lim_qtile <- qnorm(fcast_lim + (1-fcast_lim)/2)

# select heuristics
which_heu <- 'h1'

# initialize R environment
debug <- FALSE
if(debug)
	par(mfrow=c(2,4))
outfname <- sprintf("%s-%s-%d-%d-%1.1f", datafname, which_heu, max_p, max_q, errtol)
out <- file(outfname, 'w')

#- - - - - - - - END CONFIGURABLE SECTION - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Heuristics 1
# 1	Collect 50 samples
# 2	Predict l samples
# 3	If k out of l forecasts satisfy criteria, then skip sampling for k
#	times
# 4	Slide sampling window k steps forward, i.e. collect k more samples
# 5 Go back to step 3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h1 <- function() {
	success <- TRUE

	#- - - - - - - - Forecast - - - - - - - -
	fcast <- predict(best_model, n.ahead=horiz)
	z <<- c(z, fcast$pred)						# concatenate predictions to original time series
	xi <- all_data[beg_idx:(beg_idx+d-1)]		# initial values for diffinv
	z <<- diffinv(z, differences=d, xi=xi)		# recover original time series
	z_n_l <- z[(length(z)-horiz+1):length(z)]	# forecasts
	half_intrvl <- lim_qtile*fcast$se			# half of confidence interval (search 'AirPassengers' in R manual)
	lb <- z_n_l - half_intrvl					# lower bounds of confidence intervals
	ub <- z_n_l + half_intrvl					# upper bounds of confidence intervals

	#- - - - - - - - Skip sampling - - - - - - - -
	skip <- 0
	for(i in 1:horiz) {
		if(half_intrvl[i] < errtol) {
			skip <- skip + 1
			str_in <- ifelse(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i],
				' in', 'nin')
			write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f, %s %f-wide [%f, %f]",
				end_idx+i, all_data[end_idx+i], z_n_l[i], abs(all_data[end_idx+i]-z_n_l[i]),
				str_in, ub[i]-lb[i], lb[i], ub[i]), out)
		} else
			break
	}
	for(i in (skip+1):horiz) {
		str_in <- ifelse(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i],
			' in', 'nin')
		write(sprintf(" t=%3d, actual=%f, fcast=%f, err=%f, %s %f-wide [%f, %f]",
			end_idx+i, all_data[end_idx+i], z_n_l[i], abs(all_data[end_idx+i]-z_n_l[i]),
			str_in, ub[i]-lb[i], lb[i], ub[i]), out)		
	}	
	tot_skips <<- tot_skips + skip
	tot_smpls <<- end_idx

	#- - - - - - - - Resume sampling - - - - - - - -
	# sample slide_step more times before attempting to forecast again
	if(skip > 0)
		slide_step <- skip
	else
		slide_step <- horiz
	beg_idx <<- beg_idx + skip + slide_step
	if(beg_idx > max_idx)
		success <- FALSE
	z <<- all_data[beg_idx:end_idx]
	if(skip > 0) {
		z <<- c(z, z_n_l[1:skip])
		fin_data[(end_idx+1):(end_idx+skip)] <<- z_n_l[1:skip]
	}
	end_idx <<- end_idx + skip + slide_step
	if(end_idx > max_idx)
		success <- FALSE
	z <<- c(z, all_data[(end_idx-slide_step+1):end_idx])
	write(sprintf("sizeof(z)=%d", length(z)), out)
	z <<- ts(z, start=beg_idx)

	return(success)
}

#- - - - - - - - Main - - - - - - - -
tot_skips <- 0
tot_smpls <- 0

z = ts(all_data[beg_idx:end_idx], start=beg_idx)

repeat {
	if(debug) plot(z)
	
	#- - - - - - - - Power transformation - - - - - - - -
	# TODO: determine whether to perform power transformation

	#- - - - - - - - Differencing - - - - - - - -	
	d <- 2
	if(d > 0)
		z <- diff(z, differences=d)
	if(debug) plot(z)

	#- - - - - - - - Model identification and estimation - - - - - - - -
	if(debug) {
		acf(z, lag.max=50)
		pacf(z, lag.max=50)
	}

	min_aic <- 1e+6
	for(p in 0:max_p) {
		for(q in 0:max_q) {
			if(p == 0 && q == 0) next
			model <- NULL
			tryCatch(model <- arima(z, order=c(p,d,q), include.mean=TRUE),
				error=function(e) {
					if(debug) {
						writeLines(sprintf("[p=%d, q=%d] ", p, q), sep="");
						writeLines(toString(e), sep="")
					}
				})
			if(!is.null(model)) {			
				if(debug) 
					write(sprintf("[p=%d, q=%d] aic=%f, sigma=%f", p, q, model$aic, sqrt(model$sigma2)), out)
				if(model$aic < min_aic) {
					min_aic <- model$aic
					best_model <- model
					best_p <- p
					best_q <- q
				}
			}
		}
	}
	write(sprintf("Best model: p=%d, q=%d, aic=%f, sigma=%f",
		best_p, best_q, best_model$aic, sqrt(best_model$sigma2)), out)

	#- - - - - - - - Forecast and adjust next sampling window - - - - - - - -
	if(!eval(call(which_heu)))
		break
}

write(sprintf("max_p=%d, max_q=%d, tot_skips=%d, tot_smpls=%d",
	max_p, max_q, tot_skips, tot_smpls), out)

# compare raw and skip-sampled data
dev.new()
plot(ts(all_data[1:tot_smpls], start=1), main=datafname, ylab=datatype)
par(col='blue')
lines(ts(fin_data[1:tot_smpls], start=1))
par(col='black')

close(out)
