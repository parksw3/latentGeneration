
current: target
-include target.mk

##################################################################

## Defs

Sources += Makefile

ms = makestuff
Ignore += $(ms)
Makefile: $(ms)/Makefile
	touch $@
$(ms)/Makefile:
	ls ../makestuff/Makefile && /bin/ln -s ../makestuff 

-include $(ms)/os.mk

######################################################################

Sources += $(wildcard *.R)

data.Rout: data.R

mers_analysis_full_relax.Rout: mers_analysis_full_relax.R

######################################################################

### Makestuff

-include $(ms)/stepR.mk
-include $(ms)/git.mk
-include $(ms)/visual.mk
