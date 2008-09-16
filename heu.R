# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation of various heuristics
# To run this simulation, do e.g.:
# which_data=1; which_heur=0; source('heu.R')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
require('stats')
require('flexmix')
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
max_idx <- 5000
max_dsp_smpls <- 500

# select heuristics
#which_heur <- 1	# to be specified via cmd line
if(!exists("which_heur"))
	stop('Must specify which_heur')

# maximum p and q
max_p <- 4
max_q <- 4

# maximum forecasting horizon (or maximum skip samples limit in Supriyo's terminology)
mssl <- 5

# forecast limit
fcast_lim <- 0.90
lim_qtile <- qnorm(fcast_lim + (1-fcast_lim)/2)

# initialize R environment
out_fname <- sprintf("%s-%d-%d-%d-%d-%1.1f", data_fname, max_idx, which_heur, max_p, max_q, err_tol)
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
# Heuristic 1 (Supriyo's algorithm that uses forecasted values instead of interpolated values)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
init0 <- function() { #{{{
	cssl <<- 0		# CurrentSkipSampleSize
} # }}}

init1 <- function() { #{{{
	cssl <<- 0		# CurrentSkipSampleSize
} # }}}

heur0or1 <- function(selector) { #{{{
	#- - - - - - - - Step 2 - - - - - - - -
	# predict 1 sample
	tryCatch(best_model <<- ar.yw(z, order.max=1),
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
		if(selector == 0) {
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
			# Heuristics 1: keep forecasts as samples, which has already been done
		}

		fin_data[(end_idx+1):(end_idx+cssl)] <<- z[(length(z)-cssl):(length(z)-1)]
		fcast_errs <- abs(
			all_data[(end_idx+1):(end_idx+cssl)] - fin_data[(end_idx+1):(end_idx+cssl)])
		fin_errs <<- c(fin_errs, fcast_errs)

		# diagnostics
		for(i in 1:cssl) {
			write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f", 
				end_idx+i, all_data[end_idx+i], fin_data[end_idx+i], fcast_errs[i]), out)
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
	# slide sampling window one step forward
	beg_idx <<- beg_idx + horiz
	end_idx <<- end_idx + horiz
	if(end_idx > max_idx)
		return(FALSE)	
	z <<- ts(z[(length(z)-buf_len+1):length(z)], start=beg_idx)
	tot_smpls <<- end_idx
	
	# compare the original samples and our reduced samples
	#cat("\n all_data[beg_idx:end_idx]=\n   ", all_data[beg_idx:end_idx], file=out)
	#cat("\n z=", z, file=out)
	cat("\n diff=", all_data[beg_idx:end_idx]-z, "\n\n", file=out)

	return(TRUE)
} #}}}

heur0 <- function() {
	return(heur0or1(0))
}

heur1 <- function() {
	return(heur0or1(1))
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Heuristics 2
# 1	Collect 50 samples
# 2	Predict l samples
# 3	If k out of l forecasts satisfy criteria, then skip sampling for k
#	times
# 4	Slide sampling window k steps forward, i.e. collect k more samples
# 5 Go back to step 3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
init2 <- function() { #{{{
} #}}}

heur2 <- function() { #{{{
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
						writeLines(sprintf("[p=%d, q=%d] ", p, q), sep="")
						writeLines(toString(e), sep="")
					}
				})
			if(!is.null(model)) {			
				if(debug) 
					write(sprintf("[p=%d, q=%d] aic=%f, sigma=%f", p, q, model$aic, sqrt(model$sigma2)), out)
				if(model$aic < min_aic) {
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
	fcast_errs <- abs(
		all_data[(end_idx+1):(end_idx+mssl)] - z_n_l[1:mssl])

	#- - - - - - - - Skip sampling - - - - - - - -
	skip <- 0
	for(i in 1:mssl) {
		if(half_intrvl[i] < err_tol) {
			skip <- skip + 1
			fin_data[end_idx+i] <<- z_n_l[i]

			# diagnostics
			str_in <- ifelse(all_data[end_idx+i] > lb[i] && all_data[end_idx+i] < ub[i],
				' in', 'nin')
			write(sprintf("*t=%3d, actual=%f, fcast=%f, err=%f, %s %f-wide [%f, %f]",
				end_idx+i, all_data[end_idx+i], z_n_l[i], fcast_errs[i],
				str_in, ub[i]-lb[i], lb[i], ub[i]), out)
		} else
			break
	}
	# diagnostics
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

	#- - - - - - - - Resume sampling - - - - - - - -
	# sample slide_step more times before attempting to forecast again
	if(skip > 0)
		slide_step <- skip
	else
		slide_step <- mssl
	z <<- c(z[1:(buf_len+skip)], all_data[(end_idx+skip+1):(end_idx+skip+slide_step)])
	beg_idx <<- beg_idx + skip + slide_step
	end_idx <<- end_idx + skip + slide_step
	if(end_idx > max_idx)
		return(FALSE)
	z <<- ts(z[(length(z)-buf_len+1):length(z)], start=beg_idx)
	tot_smpls <<- end_idx + skip

	# compare the original samples and our reduced samples	
	#cat("\n all_data[beg_idx:end_idx]=\n       ", all_data[beg_idx:end_idx], file=out)
	#cat("\n z=", z, "\n\n", file=out)
	cat("\n diff=", all_data[beg_idx:end_idx]-z, "\n\n", file=out)

	return(TRUE)
} #}}}

#- - - - - - - - Main - - - - - - - -
fin_data <- all_data
fin_errs <- NULL

tot_skips <- 0
tot_smpls <- 0

z = ts(all_data[beg_idx:end_idx], start=beg_idx)
best_model <- NULL	# it's good to have this because
					# in case some iteration fails to find a model
					# the later iteration can use the prev model

# does not return
func_init <- sprintf("init%d", which_heur)

# returns success (boolean)
func_heur <- sprintf("heur%d", which_heur)

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

# compute Kullback-Leibler divergence
m = cbind(all_data[1:tot_smpls], fin_data[1:tot_smpls])
kl_div = KLdiv(cbind(all_data[1:tot_smpls], fin_data[1:tot_smpls]))[1,2]

# display result
write(sprintf("max_p=%d, max_q=%d, tot_skips=%d(%d), tot_smpls=%d, reductn=%f, kl_div=%g",
	max_p, max_q, tot_skips, length(fin_errs), tot_smpls, tot_skips/tot_smpls, kl_div), out)
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
if(tot_smpls <= max_dsp_smpls) {
	plot(ts(all_data[1:tot_smpls], start=1), main=data_fname, ylab=data_type)
} else {
	plot(ts(all_data[1:max_dsp_smpls], start=1), main=data_fname, ylab=data_type)
}
par(col='blue')
if(tot_smpls <= max_dsp_smpls) {
	lines(ts(fin_data[1:tot_smpls], start=1))
} else {
	lines(ts(fin_data[1:max_dsp_smpls], start=1))
}
legend("bottomright", sprintf('reduction=%f\nkl_div=%g', tot_skips/tot_smpls, kl_div))
par(col='black')
if(!debug) {
	dev.off()	# needed to flush the PDF
}

close(out)

# vim:foldmethod=marker:
