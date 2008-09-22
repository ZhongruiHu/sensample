# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation of various heuristics
# To run this simulation, do e.g.:
# which_heur='0'; which_data=1; source('heu.R')
# Additionally, max_idx, max_p and max_q can be set externally
#
# For convenience, copy and paste this to run in batch mode:
# which_heur='0'; for(which_data in 1:8) source('heu.R')
# which_heur='0a'; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=1; max_q=0; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=0; max_q=1; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=3; max_q=0; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=6; max_q=0; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=max_q=1; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=max_q=2; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=max_q=3; for(which_data in 1:8) source('heu.R')
# which_heur='1'; max_p=max_q=4; for(which_data in 1:8) source('heu.R')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
require('stats')
#require('flexmix')	# originally for KLdiv()
graphics.off()

#- - - - - - - - BEGIN CONFIGURABLE SECTION - - - - - - - -

# which data
# some data from http://www.ndbc.noaa.gov/historical_data.shtml
#which_data <- 4	# to be specified via cmd line
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

# how much data to process
buf_len <- 50
beg_idx <- 1
end_idx <- buf_len
if(!exists("max_idx"))
	max_idx <- 5000
max_dsp_smpls <- 500

# select heuristics
#which_heur <- '0'	# to be specified via cmd line
if(!exists("which_heur"))
	stop('Must specify which_heur')

# maximum p and q
if(!exists("max_p"))
	max_p <- 4
if(!exists("max_q"))
	max_q <- 4

# maximum forecasting horizon (or maximum skip samples limit in Supriyo's terminology)
mssl <- 5

# forecast limit
fcast_lim <- 0.90
lim_qtile <- qnorm(fcast_lim + (1-fcast_lim)/2)

# initialize R environment
out_fname <- sprintf("%s-%d-%s-%d-%d-%1.1f", data_fname, max_idx, which_heur, max_p, max_q, err_tol)
out <- file(out_fname, 'w')
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
#
# Heuristic 0a (Supriyo's algorithm that uses forecasted values instead of interpolated values)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
init0 <- function() { #{{{
	cssl <<- 0		# CurrentSkipSampleSize
	max_p <<- 1
	max_q <<- 0
} #}}}

init0a <- function() { #{{{
	cssl <<- 0		# CurrentSkipSampleSize
	max_p <<- 1
	max_q <<- 0	
} # }}}

heur0core <- function(selector) { #{{{
	#- - - - - - - - Step 2 - - - - - - - -
	# predict 1 sample
	tryCatch(best_model <<- ar.yw(z, order.max=max_p),
		error=function(e) {
			if(debug) {
				writeLines(sprintf("[beg_idx=%d] ", beg_idx), sep="")
				writeLines(toString(e), sep="")
			}
		})
	horiz <- cssl + 1							# past cssl forecast (past skipped samples) plus 1 current forecast
	fcast <- predict(best_model, n.ahead=horiz)	# remember cssl hasn't been incremented yet
	z <<- c(z, fcast$pred)						# concatenate predictions to original time series
	z <<- diffinv(z, differences=d, xi=z_iv)	# recover original time series
	z_n_l <- z[length(z)]						# retrieve the last forecast	
	tot_skips <<- tot_skips + cssl

	# read 1 sample
	z[length(z)] <<- all_data[end_idx+horiz]

	if(cssl > 0) {		
		if(selector == '0') {
			# Heuristics 0: use interpolated values as samples, which is done below
			#    y1 = y0 + (x1-x0)(y2-y0)/(x2-x0)
			x0 <- length(z)-horiz; y0 <- z[x0]
			x2 <- length(z); y2 <- z[x2]
			for(i in 1:cssl) {			
				x1 <- length(z)-horiz+i
				z[x1] <<- y0 + (x1-x0)*(y2-y0)/(x2-x0)
				write(sprintf(" x1=%d, y1=%f, x0=%d, y0=%f, x2=%f, y2=%f", x1, z[x1], x0, y0, x2, y2), out)
			}
		} else {
			# Heuristics 0a: keep forecasts as samples, which has already been done
		}

		fin_data[(end_idx+1):(end_idx+cssl)] <<- z[(length(z)-cssl):(length(z)-1)]

		# diagnostics
		for(i in 1:cssl) {
			err <- abs(all_data[end_idx+i] - fin_data[end_idx+i])
			write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f", 
				end_idx+i, all_data[end_idx+i], fin_data[end_idx+i], err), out)
		}
	}

	#- - - - - - - - Step 3 - - - - - - - -	
	err <- abs(z[length(z)] - z_n_l)
	if(err < err_tol) {
		cssl <<- min(cssl + 1, mssl)
		write(sprintf("@t=%3d, actual=%f, fcast=%f, err=%f, cssl=%d", 
			end_idx+horiz, z[length(z)], z_n_l, err, cssl), out)
	} else {
		cssl <<- 0
		write(sprintf(" t=%3d, actual=%f, fcast=%f, err=%f, resets cssl", 
			end_idx+horiz, z[length(z)], z_n_l, err), out)
	}

	#- - - - - - - - Step 4 - - - - - - - -
	# slide sampling window horiz steps forward
	beg_idx <<- beg_idx + horiz
	end_idx <<- end_idx + horiz
	tot_smpls <<- end_idx
	if(end_idx > max_idx)
		return(FALSE)	
	z <<- ts(z[(length(z)-buf_len+1):length(z)], start=beg_idx)
	
	# compare the original samples and our reduced samples
	#cat("\n all_data[beg_idx:end_idx]=\n   ", all_data[beg_idx:end_idx], file=out)
	#cat("\n z=", z, file=out)
	cat("\n diff=", abs(all_data[beg_idx:end_idx]-z), "\n\n", file=out)

	return(TRUE)
} #}}}

heur0 <- function() { #{{{
	return(heur0core('0'))
} #}}}

heur0a <- function() { #{{{
	return(heur0core('0a'))
} #}}}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Heuristics 1
# 1	Collect 50 samples
# 2	Forecast l samples
# 3	If k out of confidence intervals < error threshold, then skip sampling
#	for k times
# 4	Slide sampling window k steps forward, i.e. collect k more samples
# 5 Go back to step 3
#
# Heuristics 1a
# If the first forecast confidence interval > err_tol, then take that
# as the err_tol
#
# Heuristics 1b
# Combines Heuristics 1a with the using of interpolated values in place
# of forecasts
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
init1 <- function() { #{{{
} #}}}

init1a <- function() { #{{{
} #}}}

init1b <- function() { #{{{
} #}}}

heur1core <- function(selector) { #{{{
	#- - - - - - - - Model identification and estimation - - - - - - - -
	if(debug) {
		acf(z, lag.max=50)
		pacf(z, lag.max=50)
	}

	min_aic <- NULL
	best_p <- 1	# it might happen that no model can be obtained
	best_q <- 0	# in which case these are the default p and q
	for(p in 0:max_p) {
		for(q in 0:max_q) {
			if(p == 0 && q == 0) next
			model <- NULL
			tryCatch(model <- arima(z, order=c(p,d,q), include.mean=TRUE),
				error=function(e) {
					if(debug) {
						writeLines(sprintf("[p=%d, q=%d] ", p, q), sep="")
						writeLines(toString(e), sep="")
					}
				})
			if(!is.null(model)) {			
				if(debug) 
					write(sprintf("[p=%d, q=%d] aic=%f, sigma=%f", p, q, model$aic, sqrt(model$sigma2)), out)
				if(length(min_aic) == 0 || model$aic < min_aic) {
					min_aic <- model$aic
					best_model <<- model
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
	z <<- diffinv(z, differences=d, xi=z_iv)	# recover original time series
	z_n_l <- z[(length(z)-mssl+1):length(z)]	# forecasts
	half_intrvl <- lim_qtile*fcast$se			# half of confidence interval (search 'AirPassengers' in R manual)
	lb <- z_n_l - half_intrvl					# lower bounds of confidence intervals
	ub <- z_n_l + half_intrvl					# upper bounds of confidence intervals

	#- - - - - - - - Skip sampling - - - - - - - -
	# Heuristics 1a, 1b: adjust err_tol if err_tol is too small
	if(selector == '1a' || selector == '1b') {
		if(beg_idx == 1 && err_tol < half_intrvl[1]) {
			err_tol <<- half_intrvl[1]
			write(sprintf("err_tol <<- %f", err_tol), out)
		}
	}
	
	# determine how many to skip
	skip <- 0
	for(i in 1:mssl) {
		if(half_intrvl[i] < err_tol) {
			skip <- skip + 1
			fin_data[end_idx+i] <<- z_n_l[i]

			# diagnostics
			err <- abs(all_data[end_idx+i] - z_n_l[i])
			str_in <- ifelse(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i],
				' in', 'nin')
			write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f, %s %f-wide [%f, %f]",
				end_idx+i, all_data[end_idx+i], z_n_l[i], err,
				str_in, ub[i]-lb[i], lb[i], ub[i]), out)
		} else {
			break
		}
	}
	tot_skips <<- tot_skips + skip

	# diagnostics
	if(skip < mssl) {
		for(i in (skip+1):mssl) {
			err <- abs(all_data[end_idx+i] - z_n_l[i])
			str_in <- ifelse(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i],
				' in', 'nin')
			write(sprintf(" t=%3d, actual=%f, fcast=%f, err=%f, %s %f-wide [%f, %f]",
				end_idx+i, all_data[end_idx+i], z_n_l[i], err,
				str_in, ub[i]-lb[i], lb[i], ub[i]), out)
		}
	}

	#- - - - - - - - Resume sampling - - - - - - - -
	# acquire n2smpls samples before attempting to forecast again
	if(skip > 0)
		n2smpls <- skip
	else
		n2smpls <- mssl
	z <<- c(z[1:(buf_len+skip)], all_data[(end_idx+skip+1):(end_idx+skip+n2smpls)])

	# Heuristics 1b: replace forecasts with interpolated values
	if(selector == '1b') {
		if(skip > 0) {
			x0 <- buf_len; y0 <- z[x0]
			x2 <- buf_len+skip+1; y2 <- z[x2]
			for(i in 1:skip) {
				x1 <- buf_len+i
				z[x1] <<- y0 + (x1-x0)*(y2-y0)/(x2-x0)
				fin_data[end_idx+i] <<- z[x1]
			}
		}
	}
	
	# advance pointers
	beg_idx <<- beg_idx + skip + n2smpls
	end_idx <<- end_idx + skip + n2smpls
	tot_smpls <<- end_idx
	if(end_idx > max_idx)
		return(FALSE)
	z <<- ts(z[(length(z)-buf_len+1):length(z)], start=beg_idx)	

	# compare the original samples and our reduced samples	
	#cat("\n all_data[beg_idx:end_idx]=\n       ", all_data[beg_idx:end_idx], file=out)
	#cat("\n z=", z, "\n\n", file=out)
	cat("\n diff=", abs(all_data[beg_idx:end_idx]-z), "\n\n", file=out)

	return(TRUE)
} #}}}

heur1 <- function() { #{{{
	return(heur1core('1'))
} #}}}

heur1a <- function() { #{{{
	return(heur1core('1a'))
} #}}}

heur1b <- function() { #{{{
	return(heur1core('1b'))
} #}}}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Heuristics 2
# 1	Collect 50 samples
# 2	Collect and forecast l sample
# 3	If close, then skip next sample; otherwise go back to Step 2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
init2 <- function() { #{{{
} #}}}

heur2 <- function() { #{{{
	min_aic <- NULL
	best_p <- 1	# it might happen that no model can be obtained
	best_q <- 0	# in which case these are the default p and q
	for(p in 0:max_p) {
		for(q in 0:max_q) {
			if(p == 0 && q == 0) next
			model <- NULL
			tryCatch(model <- arima(z, order=c(p,d,q), include.mean=TRUE),
				error=function(e) {
					if(debug) {
						writeLines(sprintf("[p=%d, q=%d] ", p, q), sep="")
						writeLines(toString(e), sep="")
					}
				})
			if(!is.null(model)) {			
				if(debug) 
					write(sprintf("[p=%d, q=%d] aic=%f, sigma=%f", p, q, model$aic, sqrt(model$sigma2)), out)
				if(length(min_aic) == 0 || model$aic < min_aic) {
					min_aic <- model$aic
					best_model <<- model
					best_p <- p
					best_q <- q
				}
			}
		}
	}
	write(sprintf("Best model: p=%d, q=%d, aic=%f, sigma=%f",
		best_p, best_q, best_model$aic, sqrt(best_model$sigma2)), out)

	#- - - - - - - - Step 2 - - - - - - - -
	fcast <- predict(best_model, n.ahead=2)		# 1 for comparison, 1 for potential substitute for actual sample
	z <<- c(z, fcast$pred)						# concatenate predictions to original time series
	z <<- diffinv(z, differences=d, xi=z_iv)	# recover original time series
	z_n_l <- z[(length(z)-1):length(z)]			# forecasts
	z[length(z)-1] <<- all_data[end_idx+1]

	#- - - - - - - - Step 3 - - - - - - - -
	err <- abs(z[length(z)-1] - z_n_l[1])
	if(err < err_tol) {
		# diagnostics
		write(sprintf("@t=%3d, actual=%f, fcast=%f, err=%f",
			end_idx+1, z[length(z)-1], z_n_l[1], err), out)

		# update count
		tot_skips <<- tot_skips + 1	

		# take the forecast as a substitute for the actual sample
		fin_data[end_idx+2] <<- z_n_l[2]

		# diagnostics
		err <- abs(all_data[end_idx+2] - fin_data[end_idx+2])
		write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f",
			end_idx+2, all_data[end_idx+2], fin_data[end_idx+2], err), out)
	} else {	
		# diagnostics
		write(sprintf(" t=%3d, actual=%f, fcast=%f, err=%f",
			end_idx+1, z[length(z)-1], z_n_l[1], err), out)

		# update z with the actual sample
		z[length(z)] <<- all_data[end_idx+2]
	}		

	# adjust pointers
	beg_idx <<- beg_idx + 2
	end_idx <<- end_idx + 2
	tot_smpls <<- end_idx
	if(end_idx > max_idx)
		return(FALSE)
	z <<- ts(z[(length(z)-buf_len+1):length(z)], start=beg_idx)

	# compare the original samples and our reduced samples	
	#cat("\n all_data[beg_idx:end_idx]=\n       ", all_data[beg_idx:end_idx], file=out)
	#cat("\n z=", z, "\n\n", file=out)
	cat("\n diff=", abs(all_data[beg_idx:end_idx]-z), "\n\n", file=out)

	return(TRUE)
} #}}}

#- - - - - - - - Main - - - - - - - -
fin_data <- all_data

tot_skips <- 0
tot_smpls <- 0

z = ts(all_data[beg_idx:end_idx], start=beg_idx)
best_model <- NULL	# it's good to have this because
					# in case some iteration fails to find a model
					# the later iteration can use the prev model

# does not return
func_init <- sprintf("init%s", which_heur)

# returns success (boolean)
func_heur <- sprintf("heur%s", which_heur)

eval(call(func_init))

repeat {
	if(debug) plot(z)
	#cat(sprintf("\n z=", length(z)), z, "\n\n", file=out)
	
	#- - - - Power transformation - - - -
	# TODO: determine whether to perform power transformation

	#- - - - Differencing - - - -	
	d <- 2
	if(d > 0) {
		z_iv <- z[1:d]	# initial values to be used by functions later to recover the original series
		z <- diff(z, differences=d)
	}
	if(debug) plot(z)

	#- - - - Run heuristics - - - -
	if(!eval(call(func_heur)))
		break
}

# compute Kullback-Leibler divergence (wrong)
#m = cbind(all_data[1:tot_smpls], fin_data[1:tot_smpls])
#kl_div = KLdiv(cbind(all_data[1:tot_smpls], fin_data[1:tot_smpls]))[1,2]

# computer RMSE
rmse <- sqrt(mean((all_data[1:tot_smpls]-fin_data[1:tot_smpls])^2))

# display result
write(sprintf("max_p=%d, max_q=%d, tot_skips=%d, tot_smpls=%d, reductn=%f, rmse=%g",
	max_p, max_q, tot_skips, tot_smpls, tot_skips/tot_smpls, rmse), out)
cat("\n all_data=", all_data[1:tot_smpls], file=out)
cat("\n fin_data=", fin_data[1:tot_smpls], file=out)

# compare raw and skip-sampled data
if(debug) {
	dev.new()
} else {
	pdf(sprintf("%s.pdf", out_fname), width=7.5, height=7.5/3, family="Times")	
}
# must be done here after because this is a new device
# note: heavy fine-tuning
par(mar=c(2.5,2.5,1.5,0.5), mgp=c(1.5,0.5,0), cex.axis=0.8)
if(buf_len+max_dsp_smpls <= tot_smpls) {
	plot(ts(all_data[buf_len:(buf_len+max_dsp_smpls)], start=buf_len), main=data_fname, ylab=data_type)
} else {
	plot(ts(all_data[1:tot_smpls], start=1), main=data_fname, ylab=data_type)
}
par(col='blue')
if(buf_len+max_dsp_smpls <= tot_smpls) {
	lines(ts(fin_data[buf_len:(buf_len+max_dsp_smpls)], start=buf_len))
} else {
	lines(ts(fin_data[1:tot_smpls], start=1))
}
legend("bottomright", sprintf('reduction=%f\nrmse=%g', tot_skips/tot_smpls, rmse))
par(col='black')
if(!debug) {
	dev.off()	# needed to flush the PDF
}

close(out)

# vim:foldmethod=marker:
