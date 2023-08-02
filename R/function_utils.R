### function_utils.R contains general utility functions used by multiple functions/scripts:
## - get_min_age: Get minimum age from age group string, e.g. 0 from "0-9"
## - get_max_age: Get maximum age from age group string, e.g. 9 from "0-9"
## - get_reg_nb: Get region number from Google COVID-19 Open Data location_key variable , e.g. "01" from "FR_ARA_01"
## - get_reg_age_nb: Get region-age number from Google COVID-19 Open Data location_key joined to age group, e.g. "01.0-9" e.g. "FR_ARA_01.0-9"
## - get_age: Get age from region-age number from Google COVID-19 Open Data location_key joined to age group, e.g. "0-9" e.g. "FR_ARA_01.0-9"
## - get_nuts2_reg: Get NUTS-2 region code from Google COVID-19 Open Data location_key, e.g. "ARA" from "FR_ARA_01.0-9"
## - get_nuts3_reg: Get NUTS-3 region code from Google COVID-19 Open Data location_key, e.g. "FR_ARA_01" from "FR_ARA_01"
## - get_data: Download data with given url and file type
## - conv: Convolve basic serial with itself to get composite serial interval accounting for missing infection generations

## TAG COUNTRY
## Add if() section to get_reg_nb, get_reg_nb, get_age, get_nuts2_reg, and get_nuts3_reg
## Is used to extract the region number, region-age, age, NUTS-2 code, and NUTS-3 code
## from Google COVID-19 Open Data location_key.

get_min_age = function(x){
    as.numeric(sub("-.*","",sub("\\+|<","-",x)))  
} 

get_max_age = function(x){
    as.numeric(sub(".*-","",sub("\\+|<","-",x)))    
}

get_reg_nb = function(x, country){
    if (country == "FR"){
        substr(x, 8, 9)
    } else if (country == "CZ"){
        substr(x, 4, 5)
    } else if (country == "IT"){
        x
    }
}

get_reg_age_nb = function(x, country){
    if (country == "FR"){
        substr(x, 8, nchar(x))
    } else if (country == "CZ"){
        substr(x, 4, nchar(x))
    } 
}

get_age <- function(x, country){
    if (country == "FR"){
        substr(x, 11, nchar(x))
    } else if (country == "CZ"){
        substr(x, 7, nchar(x))
    }
}

get_nuts2_reg <- function(x, country){
    if (country == "FR"){
        substr(x, 4, 6)
    } else if (country == "CZ"){
        substr(x, 4, 4)
    } else if (country == "IT"){
        substr(x, 1, 4)
    }
}

get_nuts3_reg <- function(x, country){
    if (country == "FR"){
        substr(x, 1, 9)
    } else if (country == "CZ") {
        substr(x, 1, 5)
    } else if (country == "IT"){
        substr(x, 1, 5)
    }
}

get_data <- function(url,ft){
    temp <- tempfile()
    x <- download.file(url,temp)
    if (ft=="csv"){
        x <- read.csv(temp,na.strings="")
    } else if (ft=="tsv"){
        x <- read.delim(temp,na.strings="")
    }
    setDT(x)
    unlink(temp)
    return(x)
}

conv <- function(w_dens, max_days){
    # Compute log serial interval
    x <- log(w_dens)
    y <- rev(log(w_dens))
    # Compute the length of the serial interval
    n <- m <- length(x)
    
    r <- lapply(1:(n+m-1), function(k){
        i <- seq_len(n)
        i <- i[(k > (n - i))]
        w <- x[k-m+i]+y[i]
        
        total <- w[1]
        
        if (length(w)<2) return(total)
        
        for (i in 2:length(w)){
            total <- max(total, w[i]) + log(exp(total - max(total, w[i])) + 
                                                exp(w[i] - max(total, w[i])))
        }
        return(total)
    })[1 + 1:max_days] %>% unlist %>% exp
    return(r)
}

med_and_CI = function(x, l, u, f=1, d=1, method="round"){
    if (method == "signif"){
        paste0(signif(f*x, d), " (", signif(f*l, d), "-", signif(f*u, d), ")")
    } else if (method=="round"){
        paste0(round(f*x,d), " (", round(f*l,d), "-", round(f*u,d), ")")
    }
}

transp <- function(col, alpha=.5){
    res <- apply(col2rgb(col),2, function(c) 
        rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
    return(res)
}
