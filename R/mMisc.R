#' Extract RMSE and pseudo RSquared
#' 
#' Calculates Mean squared error and r-squared from the fitted and actual values. Returns each (respectively) in a 1x2 matrix.
#' 
#' @param model A model object of type \code{glm}
#' 
#' @examples
#' rmse.glm(mod.fit)
rmse.glm = function(m) {
	# m = quote(m.glm)
	
	n = nrow(m$model)
	k = ncol(m$model) - 1
	Y = m$model[,1]
	
	fitted = m$fitted.values
	epsilonhat = Y - fitted
	sigma2hat = sum(epsilonhat^2) / (n-k)
	rmse = sqrt(sigma2hat)
	rsq = cor(Y, fitted)^2
	cbind(rmse, rsq)
}


#' Quick plot of fitted values as a function of actual
#' 
#' X-Y scatterplot of fitted ~ actual with a line of equality 
#' 
#' @param model A model object of type \code{glm}
#' 
#' @examples
#' plot.avp(mod.fit)

### Another quick little function to plot fitted as a function of actuals
plot.avp = function(m) {
	x = m$model[,1]
	y = m$fitted.values
	
	greygrid2(x, y)
	points(y ~ x, type = "p", lty = 2, lwd = 1, col = pal.analytics[1])
	title(xlab = "Actual", ylab = "Fitted")
	# abline(lm(y~x), col = analytics[1])
	abline(a = 0, b = 1, lty = 2, col = "grey60")
}


#' Calculate a running summary 
#' 
#' Running summary (mean, sum, etc) of a vector over a moving window. 
#' 
#' @param vec Vector to summarize
#' @param w window (length) for moving summary
#' @param FUN summary function
#' 
#' @examples
#' # weekly moving average
#' runner(daily, w = 7, FUN = mean)
runner = function(vec, w, FUN) {
	fun = match.fun(FUN)
	end = length(vec)
	
	outvec = NULL	
	for (i in w:end) {
		outvec[(i-n+1)] = do.call(FUN, args = list(vec[(i-n+1):i]))
	}
	return(outvec)
}

#' Chomp off leading/trailing whitespace
#' 
#' Much needed function to chomp the leading/trailing whitespace from a string. 
#' 
#' just a wrapper around \code{gsub("^\\s+|\\s+$", "", <string>)}
#' 
#' @param x String to chomp()
#' 
#' @examples
#' chomp(x = "  This is an unruly string  ")
chomp = function (x) {gsub("^\\s+|\\s+$", "", x)} # desperately needed

#' True positive classification
#' 
#' Returns the rate at which true positives are identified in a classification problem given a certain threshold. Especially for logistic/poisson regression whether a model correctly classifies a binary outcome depends on the threshold. 
#' 
#' This function is primarily used to calculate values for the Receiver Operating Characteristics (ROC) curve function.
#'
#' @param threshold a value between 0 and 1 for the value at which an outcome is classified
#' @param truth a vector of [0 | 1] values representing the true classification
#' @param fitted a vector in the range [0,1] values representing the predicted classification
#' 
#' @examples
#' m0 = glm(outcome ~ predictor1 + predictor2, family = binomial(), data = train)
#' testprobs = predict(m0, newdata = test, type = "response")
#' truth = train$outcome
#' trueposrate(0.5, truth, testprobs)
#' # Most useful for understanding the rate of true positives across a range of thresholds
#' sapply(0:10 / 10, function(x) {trueposrate(x, truth = truth, fitted = testprobs)})
trueposrate = function(threshold, truth, fitted) {
	pred = fitted > threshold
        
	yes = truth == T
	truepos = truth == T & pred == T
	trueposrate = sum(truepos) / sum(yes)
	return(trueposrate)
}

#' False positive classification
#' 
#' Returns the rate at which false positives are identified in a classification problem given a certain threshold. Especially for logistic/poisson regression whether a model correctly classifies a binary outcome depends on the threshold. 
#' 
#' This function is primarily used to calculate values for the Receiver Operating Characteristics (ROC) curve function.
#'
#' @param threshold a value between 0 and 1 for the value at which an outcome is classified
#' @param truth a vector of [0 | 1] values representing the true classification
#' @param fitted a vector in the range [0,1] values representing the predicted classification
#' 
#' @examples
#' m0 = glm(outcome ~ predictor1 + predictor2, family = binomial(), data = train)
#' testprobs = predict(m0, newdata = test, type = "response")
#' truth = train$outcome
#' falseposrate(0.5, truth, testprobs)
#' # Most useful for understanding the rate of true positives across a range of thresholds
#' sapply(0:10 / 10, function(x) {falseposrate(x, truth = truth, fitted = testprobs)})
falseposrate = function(threshold, truth, fitted) {
	pred = fitted > threshold
        
	no = truth == F
	falsepos = truth == F & pred == T
	falseposrate = sum(falsepos) / sum(no)
	return(falseposrate)
}

#' Calculate x-y pairs of true/false positives for ROC charts
#' 
#' Returns a matrix of true/false positive rates which can be plotted as a ROC curve. 
#' 
#' @param truth a vector of [0 | 1] values representing the true classification
#' @param fitted a vector in the range [0,1] values representing the predicted classification
#' 
#' @examples
#' m0 = glm(outcome ~ predictor1 + predictor2, family = binomial(), data = train)
#' testprobs = predict(m0, newdata = test, type = "response")
#' roc.vals.train = roclines(train$outcome, m0$fitted)
#' roc.vals.test = roclines(train$outcome, testprobs)
roclines = function(truth, fitted, thresholds = 0:100 / 100) {
# 	threshold = 0.5
#	truth = a$match
# 	fitted = model$fitted.values

    tps = sapply(thresholds, function(x) {trueposrate(x, truth = truth, fitted = fitted)})
    fps = sapply(thresholds, function(x) {falseposrate(x, truth = truth, fitted = fitted)})
    
	cbind(fpr = rev(fps), tpr = rev(tps)) 
}


#' Convert R color strings to hex 
#' 
#' Quickly create a hex string for a color + alpha channel
#' 
#' @param basecol An R base color name
#' @param achannel (ideally) a hex value string for the alpha channel. Decimal values will be interpreted as hex
#' 
#' @export
#' @examples
#' col2hex("lightsteelblue", "70")
#' sapply(brewer.pal(n = 5, name = "BuPu"), function(x) {col2hex(x, achannel = '80')})
col2hex = function(basecol, achannel = "") {
	col2rgb(basecol) %>% as.hexmode %>% paste0(collapse = "") %>% paste0("#", ., achannel, collapse = "")
}

#' Calculate location quotients
#' 
#' Converts a data frame of counts to a data frame of location quotients.
#' 
#' @param d a data frame of arbitrary size
#' 
#' The location quotient is a proportion of proportions. The typical use case is comparing local proportions to some standard (e.g. national) proportion as a quick way of understanding how exaggerated local differences are. 
#' 
#' If the local proportion of households with children is 0.15 and the national proportion is 0.10 then the LQ is \code{0.15 / 0.10 = 1.5}. Which is to say that the local proportion is half again larger than the national proportion.
#' 
#' Will use the row/proportions as the local proportions and the column totals as the reference.
#' 
#' @export
#' @examples
#' lq_frame = lq(data)
lq = function(d) {
	D = eval(d)
	props = sweep(D, 1, apply(D, 1, sum), "/")
	marsum = apply(D, 2, sum)
	marprops = marsum / sum(marsum)

	sweep(props, 2, marprops, "/")
}
