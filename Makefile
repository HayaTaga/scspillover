SCSPILL_MODE ?= full
SCSPILL_TARGET ?= all

.PHONY: setup run clean-output

setup:
	SCSPILL_MODE=$(SCSPILL_MODE) Rscript code/00_setup.R

run:
	SCSPILL_MODE=$(SCSPILL_MODE) Rscript code/99_run_all.R $(SCSPILL_MODE) $(SCSPILL_TARGET)

clean-output:
	rm -rf output/figures/* output/tables/* output/logs/*
	mkdir -p output/figures output/tables output/tables/mc_result output/logs
	touch output/figures/.gitkeep output/tables/.gitkeep output/tables/mc_result/.gitkeep output/logs/.gitkeep
