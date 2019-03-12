library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(outbreaks)
library(rstan)

ll <- mers_korea_2015$linelist
cc <- mers_korea_2015$contacts %>% 
	group_by(to) %>%
	filter(diff_dt_onset == max(diff_dt_onset))

load("../analysis/mers_analysis_full_relax.rda")

ee <- rstan::extract(genfit)

rhat <- bayesplot::rhat(genfit)

hist(rhat)

report_data <- 
	data.frame(
		dt_report=as.numeric(ll$dt_report)-16543,
		id=ll$id
	)

Rtdata <- lapply(1:2000, function(x) {
	R <- ee$popR[x]
	phi <- ee$phi[x]
	
	
	data.frame(
		x=seq(0, 60, by=0.1),
		Rt=R * exp(- phi * seq(0, 60, by=0.1))
	)
}) %>%
	bind_rows(.id="sample") %>%
	group_by(x) %>%
	summarize(
		median=median(Rt),
		lwr=quantile(Rt, 0.025),
		upr=quantile(Rt, 0.975)
	)

gR <- ggplot(Rtdata) +
	geom_ribbon(aes(x, ymin=lwr, ymax=upr), alpha=0.2) +
	geom_line(aes(x, median)) +
	geom_vline(xintercept=quantile(ee$dt_infected[,1], 0.025), lty=2) +
	geom_vline(xintercept=quantile(ee$dt_infected[,1], 0.975), lty=2) +
	scale_x_continuous("Time (days)", expand=c(0, 0)) +
	scale_y_continuous("Reproductive number")

grep <- ggplot(report_data) +
	geom_bar(aes(dt_report), fill=NA, col="black") +
	geom_vline(xintercept=quantile(ee$dt_infected[,1], 0.025), lty=2) +
	geom_vline(xintercept=quantile(ee$dt_infected[,1], 0.975), lty=2) +
	scale_x_continuous("Time (days)", limits=c(0, 60), expand=c(0, 0)) +
	scale_y_continuous("Number of reported cases", expand=c(0, 0), limits=c(0, 25))	

g1 <- arrangeGrob(gR, grep, ncol=1)

ggsave("mers_time_series.pdf", g1, width=8, height=6)

R0 <- sapply(1:2000, function(x) {
	ee$popR[x] * exp(-ee$phi[x] * ee$dt_infected[x,1])
})

median(R0)
quantile(R0, c(0.025, 0.975))

R_onset <- sapply(1:2000, function(x) {
	ee$popR[x] * exp(-ee$phi[x] * (as.numeric(ll$dt_onset[1]) - 16543))
})

median(R_onset)
quantile(R_onset, c(0.025, 0.975))

kernel_primary <- lapply(1:2000, function(x) {
	R_i <- ee$indR[x,1]
	phi <- ee$phi[x]
	dt_infected <- ee$dt_infected[x,1]
	shape_gen <- ee$shape_gen[x]
	rate_gen <- ee$rate_gen[x]
	t <- seq(0, 30, by=0.1)
	
	data.frame(
		t=t,
		y=R_i * exp(- phi * (t + dt_infected)) * dgamma(t, shape_gen, rate_gen)
	)
}) %>%
	bind_rows(.id="sample") %>%
	group_by(t) %>%
	summarize(
		median=median(y),
		lwr=quantile(y, 0.025),
		upr=quantile(y, 0.975)
	) %>%
	mutate(
		case="Case 1"
	)

kernel_14 <- lapply(1:2000, function(x) {
	R_i <- ee$indR[x,14]
	phi <- ee$phi[x]
	dt_infected <- ee$dt_infected[x,14]
	shape_gen <- ee$shape_gen[x]
	rate_gen <- ee$rate_gen[x]
	t <- seq(0, 30, by=0.1)
	
	data.frame(
		t=t,
		y=R_i * exp(- phi * (t + dt_infected)) * dgamma(t, shape_gen, rate_gen)
	)
}) %>%
	bind_rows(.id="sample") %>%
	group_by(t) %>%
	summarize(
		median=median(y),
		lwr=quantile(y, 0.025),
		upr=quantile(y, 0.975)
	) %>%
	mutate(
		case="Case 14"
	)

kernel_avg <- lapply(1:2000, function(x) {
	R_i <- ee$indR[x,51]
	phi <- ee$phi[x]
	dt_infected <- ee$dt_infected[x,51]
	shape_gen <- ee$shape_gen[x]
	rate_gen <- ee$rate_gen[x]
	t <- seq(0, 30, by=0.1)
	
	data.frame(
		t=t,
		y=R_i * exp(- phi * (t + dt_infected)) * dgamma(t, shape_gen, rate_gen)
	)
}) %>%
	bind_rows(.id="sample") %>%
	group_by(t) %>%
	summarize(
		median=median(y),
		lwr=quantile(y, 0.025),
		upr=quantile(y, 0.975)
	) %>%
	mutate(
		case="Case 51"
	)

kernel_data <- bind_rows(kernel_primary, kernel_14, kernel_avg)

gkernel <- ggplot(kernel_data) +
	geom_ribbon(aes(t, ymin=lwr, ymax=upr), alpha=0.2) +
	geom_line(aes(t, median)) +
	scale_x_continuous("Generation (days)") +
	scale_y_continuous("Kernel density") +
	facet_wrap(~case) +
	theme(
		strip.background = element_blank()
	)

ggsave("kernel.pdf", gkernel, width=8, height=3)

gdistribution <- lapply(1:2000, function(x){
	shape_gen <- ee$shape_gen[x]
	rate_gen <- ee$rate_gen[x]
	t <- seq(0, 30, by=0.1)
	
	data.frame(
		t=t,
		y=dgamma(t, shape_gen, rate_gen)
	)
}) %>%
	bind_rows(.id="sample") %>%
	group_by(t) %>%
	summarize(
		median=median(y),
		lwr=quantile(y, 0.025),
		upr=quantile(y, 0.975)
	) %>%
	mutate(
		type="Generation interval"
	)

idistribution <- lapply(1:2000, function(x){
	shape_inf <- ee$shape_inf[x]
	rate_inf <- ee$rate_inf[x]
	t <- seq(0, 30, by=0.1)
	
	data.frame(
		t=t,
		y=dgamma(t, shape_inf, rate_inf)
	)
}) %>%
	bind_rows(.id="sample") %>%
	group_by(t) %>%
	summarize(
		median=median(y),
		lwr=quantile(y, 0.025),
		upr=quantile(y, 0.975)
	) %>%
	mutate(
		type="Infectious period"
	)

rdistribution <- lapply(1:2000, function(x){
	rate_report <- ee$rate_report[x]
	t <- seq(0, 30, by=0.1)
	
	data.frame(
		t=t,
		y=dexp(t, rate_report)
	)
}) %>%
	bind_rows(.id="sample") %>%
	group_by(t) %>%
	summarize(
		median=median(y),
		lwr=quantile(y, 0.025),
		upr=quantile(y, 0.975)
	) %>%
	mutate(
		type="Reporting period"
	)

ddata <- bind_rows(gdistribution, idistribution, rdistribution)

gdistr <- ggplot(ddata) +
	geom_ribbon(aes(t, ymin=lwr, ymax=upr), alpha=0.2) +
	geom_line(aes(t, median)) +
	scale_x_continuous("Time (days)") +
	scale_y_continuous("Density") +
	facet_wrap(~type) +
	theme(
		strip.background = element_blank()
	)

ggsave("distribution.pdf", gdistr, width=8, height=3)

ggen1 <- (gdistr %+% filter(ddata, type=="Generation interval")) +
	geom_histogram(data=cc, aes(diff_dt_onset, y=..density..), bins=30, fill=NA, col="black")

ggsave("generation.pdf", ggen1, width=6, height=4)
