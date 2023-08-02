## Function to create the sts object used for hhh4 models
create_sts <- function(country, counts, pop, map, C = NULL,
                       by = c("districts", "agegroups", "all", "none"),
                       agegroups = rep(1,10),
                       timeRange = c("2020-05-13", "2022-02-10"),
                       flatten = FALSE)
{
    by <- match.arg(by)
    ## Re-organise the rows in map
    regs <- get_reg_nb(colnames(counts), country)
    map <- map[match(regs, map$reg_nb),]
    ## Overwrite map row names with region names
    row.names(map) <- colnames(counts)[match(map$reg_nb,regs)]
    ## Make sure map is SpatialPolygonsDataFrame
    if (class(map)[1] != "SpatialPolygonsDataFrame"){
        map <- as(map,"Spatial")    
    }
    ## subset counts to time range
    stopifnot(!is.na(
        timeRangeIdx <- match(timeRange, dimnames(counts)[[1L]])
    ))
    counts <- counts[timeRangeIdx[1L]:timeRangeIdx[2L], , , drop = FALSE]
    dates <- as.Date(dimnames(counts)[[1L]])
    year_start <- min(lubridate::year(dates))
    start <- c(year_start, min(dates) - as.Date(paste0(year_start, "-01-01")))
    
    ## aggregate groups
    if (by %in% c("agegroups", "all") && !is.null(agegroups)) {
        counts <- aggregateCountsArray(counts = counts, dim = 3, grouping = agegroups)
        pop <- aggregateCountsArray(counts = pop, dim = 2, grouping = agegroups)
    }
    
    ## neighbourhood structure
    if (by %in% c("districts", "all")) {
            neighbourhood <- nbOrder(poly2adjmat(map), maxlag = 12)
            neighbourhood[neighbourhood == 0] <- 50
            diag(neighbourhood) <- 0
    } else if (by == "agegroups") {
        neighbourhood <- my.contactmatrix(C, colSums(pop), grouping = agegroups)
    } else { # no stratification
        neighbourhood <- matrix(NA, 1L, 1L)
    }
    ## constructor function
    makests <- function (observed, populationFrac, neighbourhood, ...) {
        new("sts", start = start, freq = 365, observed = observed, 
            populationFrac = populationFrac, neighbourhood = neighbourhood, 
            epoch = as.numeric(dates), ...)
    }
    
    ## create "sts" object for the chosen level of stratification
    if (by == "none") {
        observed <- cbind(rowSums(counts))
        colnames(observed) <- substr(dimnames(counts)[[2]][1],1,2)
        makests(
            observed = observed,
            populationFrac = 1, neighbourhood = neighbourhood
        )
    } else if (by == "districts") {
        makests(
            observed = rowSums(counts, dims = 2L),
            populationFrac = pop, map = map, 
            neighbourhood = neighbourhood
        )
    } else if (by == "agegroups") {
        makests(
            observed = apply(counts, c(1,3), sum),
            populationFrac = prop.table(colSums(pop)),
            neighbourhood = neighbourhood
        )
    } else if (flatten) { # where districts vary faster than age groups
        ngroups <- dim(counts)[[3L]]
        nregions <- dim(counts)[[2L]]
        ## replicate 'neighbourhood' block ngroup times
        ## replicate columns
        neighbourhood <- rep.int(neighbourhood, ngroups)
        dim(neighbourhood) <- c(nregions, nregions * ngroups)
        ## replicate rows
        neighbourhood <- do.call("rbind", rep_len(
            list(as.name("neighbourhood")), ngroups))
        makests(
            observed = as.data.frame(counts),
            populationFrac = c(pop),
            neighbourhood = neighbourhood
        )
    } else { # list of by="districts" sts objects by age group
        sapply(dimnames(counts)[[3L]], function (g) {
            map$POPULATION <- pop[, g]
            makests(
                observed = counts[, , g],
                populationFrac = prop.table(pop[, g]),
                map = map, neighbourhood = neighbourhood
            )
        }, simplify = FALSE, USE.NAMES = TRUE)
    }
    
}

my.contactmatrix <- function(C, pop, grouping = NULL, normalize = FALSE) 
{
    if (!is.null(grouping)) {
        weights <- pop/sum(pop)
        C <- aggregateC(C = C, grouping = grouping, weights = weights)
        attr(C, "agedistri") <- 
            aggregateCountsArray(t(weights), 
                                 dim = 2, grouping = grouping, sort = TRUE)[1L, ]
    }
    if (normalize) 
        C/rowSums(C)
    else C
}
