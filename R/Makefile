R_SRC_DIR  = $(addprefix , R)
ifdef FILE
R_FILE = $(R_SRC_DIR)/$(FILE)
else
R_FILE ?= $(foreach sdir, $(R_SRC_DIR), $(wildcard $(sdir)/*.R))
endif

LOG_FILE := $(patsubst %.R, %.log, $(R_FILE))

$(info ************************************)
$(info R Source directory:     $(R_SRC_DIR))
$(info R Source files:         $(R_FILE))
$(info R Log files:            $(LOG_FILE))
$(info ************************************)

.PHONY: all

all: R/synthetic_data.log R/site_replacement.log 

R/synthetic_data.log: R/synthetic_data.R
	@echo "Compiling the analyses"
	$(info Source = $<; Destination = $@)
	cd R; Rscript $(notdir $<) > $(notdir $@)

R/site_replacement.log: R/site_replacement.R R/synthetic_data.log
	@echo "Compiling the analyses"
	$(info Source = $<; Destination = $@)
	cd R; Rscript $(notdir $<) > $(notdir $@)

# $(LOG_FILE) : %.log : %.R
# 	@echo "Compiling the analyses"
# 	$(info Source = $<; Destination = $@)
# 	cd R; Rscript $(notdir $<) > $(notdir $@)

# all: $(R_FILE)
# 	 @echo "Compiling analyses"
# 	 cd R; Rscript $(notdir $<)
