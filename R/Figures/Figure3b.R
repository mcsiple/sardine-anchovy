# Figure 3b
load(here::here("R/DataCleanup/RAM_Barange_FAO_States.RData")) # RBF (missing years filled w MARSS states)

# Distributions for null vs. observed data
# Load data
obs <- get_obs(dat = RBF,
               dsource = "Barange",
               reg = "California",
               var = "landings")
obstest <- get_wmr(obs$std_anchovy, 
                   obs$std_sardine)

xx <- get_surrogates(obs = obs,nsurrogates = 50) 
nulltest1 = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,1],
                    sardine.ts=xx$Sardine.surrogates[,1])
nulltest2 = get_wmr(anchovy.ts=xx$Anchovy.surrogates[,2],
                    sardine.ts=xx$Sardine.surrogates[,2])


# Functions used here are in the NullModel.R file
giantnull <- get_large_null(dat = RBF,
                            dsource = "Barange",
                            reg = "California",
                            var = "landings",
                            nsims=50)
str(giantnull)

distfig <- function(null, observed, return.dataframe = F){
  #' @description plots a histogram of observed WMRs vs. "null" WMR
  #' @param null is a list of three WMR matrices, each of which represents one time window (from a null surrogate)
  #' @param obs is a list of three WMR matrices, each of which represents one time window (observed time series)
  
  nl <- data.frame(ID = rep(names(null), sapply(null, length)),
                   wmrdens = unlist(null))
  obsl <- data.frame(ID = rep(names(observed), sapply(observed, length)),
                     wmrdens = unlist(observed))
  nl$ID_ord = factor(nl$ID, levels=c('less.than.5','five.ten','ten.plus'))
  obsl$ID_ord = factor(obsl$ID, levels=c('less.than.5','five.ten','ten.plus'))
  
  nl$ID_ord <- recode(nl$ID_ord, less.than.5 = "< 5 yr", five.ten = "5-10 yr",ten.plus="10+ yr" )
  obsl$ID_ord <- recode(obsl$ID_ord, less.than.5 = "< 5 yr", five.ten = "5-10 yr",ten.plus="10+ yr" )
  
  ggplot(nl,aes(x=wmrdens)) + 
    geom_density(fill="grey",colour="darkgrey",alpha=0.5) + 
    geom_density(data=obsl,fill="#41b6c4",colour="#41b6c4",alpha=0.5) +
    facet_wrap(~ID_ord,ncol = 1) +
    ylab("Density") +
    xlab("Wavelet modulus ratio (WMR)") +
    theme_classic(base_size=14)
  
  if(return.dataframe == TRUE){
    nl$nv <- "null"
    obsl$nv <- "obs"
    all <- dplyr::bind_rows(nl,obsl)
    return(all)
  }
}

test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)

# Need to get data into the following format for analysis:
# Year Sardine.est Anchovy.est      region datasource variable
newmelt <- alldat %>%  
  select(datasource,stock,year,ssb,rec,landings,sp, region, subregion) %>%  
  melt(id.vars = c("datasource","stock","region", "subregion","year"))


# Get a giant dataframe of all the WMRs for each of the regions.
stocks <- alldat %>% 
  distinct(sp, stock, region, datasource) # get unique stock-region combos
n <- nrow(stocks) # there should be 15 combos of region and datasource

compares <- alldat %>% 
  distinct(region, datasource)


# Landings ----------------------------------------------------------------
landings.wmrs <- vector()
for(s in 1:nrow(compares)){ 
  tryCatch ({
    obs <- get_obs(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "landings")
    obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
    giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "landings",nsims=50)
    null=giantnull
    observed=obstest
    test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)
    test$region <- compares$region[s]
    test$datasource <- compares$datasource[s]
    test$variable <- "landings"
    landings.wmrs <- rbind(landings.wmrs,test)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# Biomass -----------------------------------------------------------------
ssb.wmrs <- vector()
for(s in 1:nrow(compares)){ 
  tryCatch ({
    obs <- get_obs(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "ssb")
    obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
    giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "ssb",nsims=50)
    null=giantnull
    observed=obstest
    test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)
    test$region <- compares$region[s]
    test$datasource <- compares$datasource[s]
    test$variable <- "ssb"
    ssb.wmrs <- rbind(ssb.wmrs,test)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# Rec -----------------------------------------------------------------
rec.wmrs <- vector()
for(s in 1:nrow(compares)){ 
  tryCatch ({
    obs <- get_obs(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "rec")
    obstest <- get_wmr(obs$std_anchovy, obs$std_sardine)
    giantnull <- get_large_null(dat = RBF,dsource = compares$datasource[s],reg = compares$region[s],var = "rec",nsims=50)
    null=giantnull
    observed=obstest
    test <- distfig(null = giantnull,observed = obstest,return.dataframe = T)
    test$region <- compares$region[s]
    test$datasource <- compares$datasource[s]
    test$variable <- "rec"
    rec.wmrs <- rbind(rec.wmrs,test)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}





# Figure 3b_top ---------------------------------------------------------------
landings.wmrs %>% subset(datasource=="Barange" & nv=="obs") %>% 
  ggplot(aes(x=wmrdens,colour=region,fill=region)) + 
  geom_density(alpha=0.5,lwd=1.2,trim=F) + 
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  facet_wrap(~ID_ord,nrow=1,scales = "free_y") +
  ylab("Density") +
  xlab("Wavelet modulus ratio (WMR)") +
  theme_classic(base_size=14) %+replace% theme(strip.background  = element_blank()) +
  ggtitle("Landings")


# Figure 3b_bottom ---------------------------------------------------------------
ssb.wmrs %>% subset(datasource=="Barange" & nv=="obs") %>% 
  ggplot(aes(x=wmrdens,colour=region,fill=region)) + 
  geom_density(alpha=0.5,lwd=1.2,trim=F) + 
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  facet_wrap(~ID_ord,nrow=1,scales = "free_y") +
  ylab("Density") +
  xlab("Wavelet modulus ratio (WMR)") +
  scale_y_reverse() +
  theme_classic(base_size=14) %+replace% theme(strip.background  = element_blank()) +
  ggtitle("Biomass")


# These two plots were combined and edited in Adobe CS