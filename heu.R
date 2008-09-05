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

beg_idx <- 1
end_idx <- 50
max_idx <- 300

# maximum p and q
max_p <- 4
max_q <- 4

# maximum forecasting horizon (or maximum skip samples limit in Supriyo's terminology)
mssl <- 5

# forecast limit
fcast_lim <- 0.90
lim_qtile <- qnorm(fcast_lim + (1-fcast_lim)/2)

# select heuristics
which_heur <- '0'

# initialize R environment
outfname <- sprintf("%s-%d-%s-%d-%d-%1.1f", datafname, max_idx, which_heur, max_p, max_q, errtol)
out <- file(outfname, 'w')
debug <- FALSE
if(debug) {
	par(mfrow=c(2,4))
}

#- - - - - - - - END CONFIGURABLE SECTION - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Heuristics 0 (Supriyo's algorithm)
# 1	Collect r samples, and set
#	CurrentSkipSampleSize = SkipSamples = 0,
# 2	Read and predict next sample
# 3	If actual sample is close to prediction, set
#	++CurrentSkipSampleSize (but capped at MaximumSkipSamplesLimit),
#	SkipSamples = CurrentSkipSampleSize,
#	and go to step 4;
#	else set
#	CurrentSkipSampleSize = SkipSamples = 0,
#	and go to step 2
# 4	Repeat {skip 1 sample, --SkipSamples} until SkipSamples == 0
# 5 Go to step 2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
init0 <- function() { #{{{
	cssl <<- 0		# CurrentSkipSampleSize
} # }}}

heur0 <- function() { #{{{
	#- - - - - - - - Step 2 - - - - - - - -
	# predict 1 sample
	best_model <- ar.yw(z, order.max=1)
	horiz <- cssl + 1							# past cssl forecast (skipped samples), 1 current forecast
	fcast <- predict(best_model, n.ahead=horiz)	# remember cssl hasn't been incremented yet
	z <<- c(z, fcast$pred)						# concatenate predictions to original time series
	xi <- all_data[beg_idx:(beg_idx+d-1)]		# initial values for diffinv
	z <<- diffinv(z, differences=d, xi=xi)		# recover original time series
	z_n_l <- z[(length(z)-horiz+1):length(z)]	# retrieve the forecasts
	if(cssl > 0) {
		# store forecasts
		fin_data[(end_idx+1):(end_idx+cssl)] <<- z_n_l[1:cssl]
		fcast_errs <- abs(all_data[(end_idx+1):(end_idx+cssl)] - z_n_l[1:cssl])
		fin_errs <<- c(fin_errs, fcast_errs)

		# diagnostics
		for(i in 1:cssl) {
			write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f", 
				end_idx+i, all_data[end_idx+i], z_n_l[i], fcast_errs[i]), out)	
		}
	}

	# read 1 sample
	beg_idx <<- beg_idx + cssl + 1				# slide sampling window one step forward
	end_idx <<- end_idx + cssl + 1
	if(end_idx > max_idx)
		return(FALSE)
	z <<- ts(all_data[beg_idx:end_idx], start=beg_idx)

	#- - - - - - - - Step 3 - - - - - - - -	
	err <- abs(z[length(z)] - z_n_l[length(z_n_l)])
	if(err < errtol) {
		cssl <<- min(cssl + 1, mssl)
		write(sprintf("@t=%3d, actual=%f, fcast=%f, err=%f, cssl=%d", 
			end_idx, z[length(z)], z_n_l[length(z_n_l)], err, cssl), out)
	} else {		
		cssl <<- 0
		write(sprintf(" t=%3d, actual=%f, fcast=%f, err=%f, reseting cssl", 
			end_idx, z[length(z)], z_n_l[length(z_n_l)], err), out)
		return(TRUE)
	}

	#- - - - - - - - Step 4 - - - - - - - -
	tot_skips <<- tot_skips + cssl
	tot_smpls <<- end_idx + cssl

	return(TRUE)
} #}}}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Heuristics 1
# 1	Collect 50 samples
# 2	Predict l samples
# 3	If k out of l forecasts satisfy criteria, then skip sampling for k
#	times
# 4	Slide sampling window k steps forward, i.e. collect k more samples
# 5 Go back to step 3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
init1 <- function() { #{{{
} #}}}

heur1 <- function() { #{{{
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

	#- - - - - - - - Forecast - - - - - - - -
	fcast <- predict(best_model, n.ahead=mssl)
	z <<- c(z, fcast$pred)						# concatenate predictions to original time series
	xi <- all_data[beg_idx:(beg_idx+d-1)]		# initial values for diffinv
	z <<- diffinv(z, differences=d, xi=xi)		# recover original time series
	z_n_l <- z[(length(z)-mssl+1):length(z)]	# forecasts
	half_intrvl <- lim_qtile*fcast$se			# half of confidence interval (search 'AirPassengers' in R manual)
	lb <- z_n_l - half_intrvl					# lower bounds of confidence intervals
	ub <- z_n_l + half_intrvl					# upper bounds of confidence intervals
	fcast_errs <- abs(
		all_data[(end_idx+1):(end_idx+mssl)] - z_n_l[1:mssl])

	#- - - - - - - - Skip sampling - - - - - - - -
	skip <- 0
	for(i in 1:mssl) {
		if(half_intrvl[i] < errtol) {
			skip <- skip + 1
			str_in <- ifelse(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i],
				' in', 'nin')
			write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f, %s %f-wide [%f, %f]",
				end_idx+i, all_data[end_idx+i], z_n_l[i], fcast_errs[i],
				str_in, ub[i]-lb[i], lb[i], ub[i]), out)
		} else
			break
	}
	for(i in (skip+1):mssl) {
		str_in <- ifelse(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i],
			' in', 'nin')
		write(sprintf(" t=%3d, actual=%f, fcast=%f, err=%f, %s %f-wide [%f, %f]",
			end_idx+i, all_data[end_idx+i], z_n_l[i], fcast_errs[i],
			str_in, ub[i]-lb[i], lb[i], ub[i]), out)
	}
	if(skip > 0)
		fin_errs <<- c(fin_errs, fcast_errs[1:skip])
	tot_skips <<- tot_skips + skip
	tot_smpls <<- end_idx + skip

	#- - - - - - - - Resume sampling - - - - - - - -
	# sample slide_step more times before attempting to forecast again
	if(skip > 0)
		slide_step <- skip
	else
		slide_step <- mssl
	beg_idx <<- beg_idx + skip + slide_step
	z <<- all_data[beg_idx:end_idx]
	if(skip > 0) {
		z <<- c(z, z_n_l[1:skip])
		fin_data[(end_idx+1):(end_idx+skip)] <<- z_n_l[1:skip]
	}
	end_idx <<- end_idx + skip + slide_step
	if(end_idx + mssl > max_idx)
		return(FALSE)
	z <<- c(z, all_data[(end_idx-slide_step+1):end_idx])
	write(sprintf("sizeof(z)=%d", length(z)), out)
	z <<- ts(z, start=beg_idx)

	return(TRUE)
} #}}}

#- - - - - - - - Main - - - - - - - -
fin_data <- all_data
fin_errs <- NULL

tot_skips <- 0
tot_smpls <- 0

z = ts(all_data[beg_idx:end_idx], start=beg_idx)

# does not return
func_init <- sprintf("init%s", which_heur)

# returns success (boolean)
func_heur <- sprintf("heur%s", which_heur)

eval(call(func_init))

repeat {
	if(debug) plot(z)
	
	#- - - - Power transformation - - - -
	# TODO: determine whether to perform power transformation

	#- - - - Differencing - - - -	
	d <- 2
	if(d > 0)
		z <- diff(z, differences=d)
	if(debug) plot(z)

	#- - - - Run heuristics - - - -
	if(!eval(call(func_heur)))
		break
}

# display result
write(sprintf("max_p=%d, max_q=%d, tot_skips=%d, tot_smpls=%d, reductn=%f, max_err=%f(%d)",
	max_p, max_q, tot_skips, tot_smpls, tot_skips/tot_smpls, max(fin_errs), length(fin_errs)), out)

# compare raw and skip-sampled data
if(debug) {
	dev.new()
} else {
	pdf(sprintf("%s.pdf", outfname))
}
plot(ts(all_data[1:tot_smpls], start=1), main=datafname, ylab=datatype)
par(col='blue')
lines(ts(fin_data[1:tot_smpls], start=1))
par(col='black')
if(!debug) {
	dev.off()	# needed to flush the PDF
}

close(out)

# vim:foldmethod=marker:
