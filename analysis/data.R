library(dplyr)
library(outbreaks)

ll <- mers_korea_2015$linelist
d0 <- as.numeric(min(ll$dt_start_exp, na.rm=TRUE))

## Who has multiple contacts reported?
print(mers_korea_2015$contacts %>% 
	group_by(to) %>% 
	mutate(n=n()) %>% 
	ungroup %>% 
	filter(n > 1)
)

## For now, just keep the earliest contact for these few

cc <- (mers_korea_2015$contacts
	%>% group_by(to)
	%>% filter(diff_dt_onset == max(diff_dt_onset))
)
