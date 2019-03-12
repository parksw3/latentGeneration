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
flevel <- ll$id

standata <- data_frame(
	patient=as.numeric(factor(ll$id, level=flevel)),
	from=as.numeric(factor(cc$from[match(flevel, cc$to)], level=flevel)),
	dt_start_exp=as.numeric(ll$dt_start_exp)-16543,
	dt_end_exp=as.numeric(ll$dt_end_exp)-16543,
	dt_onset=as.numeric(ll$dt_onset)-16543,
	dt_report=as.numeric(ll$dt_report)-16543
) %>%
	mutate(
		dt_start_exp=ifelse((dt_start_exp <= dt_onset[from] & !is.na(from)) | is.na(dt_start_exp), dt_onset[from], dt_start_exp),
		dt_end_exp=pmax(pmin(dt_end_exp, dt_onset, na.rm=TRUE), dt_start_exp, na.rm=TRUE),
		dt_start_exp=dt_start_exp + ifelse(dt_start_exp==dt_end_exp, -1, 0),
		dt_start_exp=ifelse(is.na(dt_start_exp), 23, dt_start_exp), 
		from=ifelse(is.na(from), -1, from),
		missing_onset=as.numeric(is.na(dt_onset)),
		missing_end_exp=as.numeric(is.na(dt_end_exp)),
		dt_onset=ifelse(is.na(dt_onset), -1, dt_onset),
		dt_end_exp=ifelse(is.na(dt_end_exp), -1, dt_end_exp)
	) %>%
	arrange(missing_onset, missing_end_exp) %>%
	mutate(
		missing_onset=cumsum(missing_onset),
		missing_end_exp=cumsum(missing_end_exp)
	) %>%
	as.list %>%
	c(npatient=length(flevel), 
	  n_missing_onset=sum(.$missing_onset > 0), 
	  n_missing_end_exp=sum(.$missing_end_exp > 0),
	  prior_only=0)

rt <- stanc(file="../model/full_model_relax.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

genfit <- sampling(sm, standata, chain=4, iter=2000, seed=101)

save("genfit", file="mers_analysis_full_relax.rda")
