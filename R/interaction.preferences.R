#' @title Estimate species' interaction preferences based on an interaction matrix and their relative abundances
#' @param adj The adjacency matrix
#' @param row.abun The abundances (or relative abundances) of row species
#' @param col.abun The abundances (or relative abundances) of column species
#' @param include.zeros Flag to decide whether or not unobserved interactions should have their preferences estimated
#' @export
`interaction.preferences` <- function(adj, row.abun=NULL, col.abun=NULL, include.zeros=FALSE){

	#################
	# let's run a few defensive checks for quality control
	#################

	# does the adjacency matrix have row/column names?
	if(is.null(rownames(adj)))
		rownames(adj) <- as.character(1:nrow(adj))
	if(is.null(colnames(adj)))
		colnames(adj) <- as.character(1:ncol(adj))

	# what about the row abundances?
	if(!is.null(row.abun)){
		# do the sizes of everything match up? if not, this is a deal breaker
		if(nrow(adj) != length(row.abun))
			stop('Your adjacency matrix and row abundances are of different size.')
		# check whether there are names for the abundances since we use them below
		if(is.null(names(row.abun))){
			warning('Your row abundances are not named. Will assume they follow the matrix ordering.')
			names(row.abun) <- rownames(adj)
		}
	}
	# and the column abundances?
	if(!is.null(col.abun)){
		# do the sizes of everything match up? if not, this is a deal breaker
		if(ncol(adj) != length(col.abun))
			stop('Your adjacency matrix and column abundances are of different size.')
		# check whether there are names for the abundances since we use them below
		if(is.null(names(col.abun))){
			warning('Your column abundances are not named. Will assume they follow the matrix ordering.')
			names(col.abun) <- colnames(adj)
		}
	}

	#################
	# now that that's finished, it's peanut butter jelly time!
	#################

	# select the data of interest depending on whether or not we wish to consider zero-frequency interactions as informative
	if(include.zeros){
		frequencies <- adj[which(!is.na(adj))]
		row.species <- rownames(adj)[which(!is.na(adj),arr.ind=T)[,"row"]]
		col.species <- colnames(adj)[which(!is.na(adj),arr.ind=T)[,"col"]]
	}else{
		frequencies <- adj[which(adj>0)]
		row.species <- rownames(adj)[which(adj>0,arr.ind=T)[,"row"]]
		col.species <- colnames(adj)[which(adj>0,arr.ind=T)[,"col"]]
	}

	# turn the above into a data frame
	adj.data <- data.frame(frequency = frequencies, row.species = row.species, col.species = col.species)

	# add the abudances if desired
	if(!is.null(row.abun)) adj.data$row.abun <- row.abun[row.species]
	if(!is.null(col.abun)) adj.data$col.abun <- col.abun[col.species]

	# create a container for all this madness
	obj <- list()

	# estimate the underlying model
	if(is.null(row.abun) & is.null(col.abun)){
		obj$model <- glm(frequency~0 + row.species + col.species, family=poisson, data=adj.data)
		obj$row.abun.est <- exp(obj$model$coefficients[paste0("row.species",rownames(adj))])
		obj$row.abun.est <- obj$row.abun.est/sum(obj$row.abun.est)
		names(obj$row.abun.est) <- rownames(adj)

		obj$model <- glm(frequency~0 + col.species + row.species, family=poisson, data=adj.data)
		obj$col.abun.est <- exp(obj$model$coefficients[paste0("col.species",colnames(adj))])
		obj$col.abun.est <- obj$col.abun.est/sum(obj$col.abun.est)
		names(obj$col.abun.est) <- colnames(adj)
	}else{
		if(is.null(row.abun)){
			obj$model <- glm(frequency~0 + row.species + offset(log(col.abun)), family=poisson, data=adj.data)

			obj$row.abun.est <- exp(obj$model$coefficients[paste0("row.species",rownames(adj))])
			obj$row.abun.est <- obj$row.abun.est/sum(obj$row.abun.est)
			names(obj$row.abun.est) <- rownames(adj)
		}else{
			if(is.null(col.abun)){
				obj$model <- glm(frequency~0 + col.species + offset(log(row.abun)), family=poisson, data=adj.data)

				obj$col.abun.est <- exp(obj$model$coefficients[paste0("col.species",colnames(adj))])
				obj$col.abun.est <- obj$col.abun.est/sum(obj$col.abun.est)
				names(obj$col.abun.est) <- colnames(adj)
			}else{
				obj$model <- glm(frequency~1 + offset(log(row.abun)+log(col.abun)), family=poisson, data=adj.data)
			}
		}
	}

	# the preferences are actually the residuals of the model above
	# but we want to return them in matrix form just like the input
	obj$log.gamma <- matrix(NA,nrow(adj),ncol(adj))
	rownames(obj$log.gamma) <- rownames(adj)
	colnames(obj$log.gamma) <- colnames(adj)
	if(include.zeros){
		obj$log.gamma[which(!is.na(adj))] <- resid(obj$model)
	}else{
		obj$log.gamma[which(adj>0)] <- resid(obj$model)
	}

	return(obj)
}

	
# data$host <- as.factor(which(!is.na(net),arr.ind=T)[,"row"])
# data$host.abun <- h.abun[which(!is.na(net),arr.ind=T)[,"row"]]
# data$para <- as.factor(which(!is.na(net),arr.ind=T)[,"col"])
# data$para.abun <- p.abun[which(!is.na(net),arr.ind=T)[,"col"]]

##############################################################################
##############################################################################
### peanut butter jelly time!
##############################################################################
##############################################################################
#
# estimate the interaction "preferences" depending on various assumptions
#
# note that the following requires that the interaction data are integers (and not 0/1)
#
##############################################################################

# # estimate the preferences *and* best fit abundances according to the observed interactions
# model1 <- glm(frequency~1 + host + para, family=poisson, data=data)
# prefs1 <- resid(model1)

# # estimate the preferences *and* best fit host abundances according to the observed interactions and observed parasite abundances
# model2 <- glm(frequency~1 + host + offset(log(para.abun)), family=poisson, data=data)
# prefs2 <- resid(model2)

# # estimate the preferences *and* best fit parasite abundances according to the observed interactions and observed host abundances
# model3 <- glm(frequency~1 + para + offset(log(host.abun)), family=poisson, data=data)
# prefs3 <- resid(model3)

# # estimate the preferences according to the observed interactions, observed host abundances, and observed parasite abundances
# model4 <- glm(frequency~1+offset(log(host.abun)+log(para.abun)), family=poisson, data=data)
# prefs4 <- resid(model4)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
