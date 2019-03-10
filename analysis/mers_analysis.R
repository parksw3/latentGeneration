library(dplyr)
library(outbreaks)
library(rstan)

ll <- mers_korea_2015$linelist
cc <- mers_korea_2015$contacts %>% 
	group_by(to) %>%
	filter(diff_dt_onset == max(diff_dt_onset))

mers_korea_2015$contacts %>% 
	group_by(to) %>% 
	mutate(n=length(to)) %>% 
	ungroup %>% 
	filter(n > 1)

## factor level
known <- unique(c(cc$from, cc$to))

standata <- data_frame(
	patient=as.numeric(factor(known, level=known)),
	from=as.numeric(factor(cc$from[match(known, cc$to)], level=known)),
	dt_start_exp=as.numeric(ll$dt_start_exp[match(known, ll$id)]) - 16543,
	dt_end_exp=as.numeric(ll$dt_end_exp[match(known, ll$id)]) - 16543,
	dt_onset=as.numeric(ll$dt_onset[match(known, ll$id)]) - 16543,
	primary=as.numeric(known=="SK_1")
) %>%
	mutate(
		dt_start_exp=ifelse((dt_start_exp <= dt_onset[from] & !is.na(from)) | is.na(dt_start_exp), dt_onset[from], dt_start_exp),
		dt_end_exp=pmax(pmin(dt_end_exp, dt_onset, na.rm=TRUE), dt_start_exp, na.rm=TRUE),
		dt_start_exp=dt_start_exp + ifelse(dt_start_exp==dt_end_exp, -1, 0),
		from=ifelse(is.na(from), 0, from)
	) %>%
	as.list %>%
	c(t_censor=60, npatient=length(known), prior_only=0)

rt <- stanc(file="../model/latent_model.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

genfit <- sampling(sm, standata, chain=4, iter=4000, seed=101,
				   control=list(max_treedepth=12, adapt_delta=0.9))

save("standata", "genfit", file="mers_analysis.rda")
