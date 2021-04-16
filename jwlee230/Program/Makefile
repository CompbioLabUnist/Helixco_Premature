# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = premature:latest
CPUS = 10
MEMS = 300G
PWD := $(shell pwd)
TOOLS = 
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data --volume $(abspath Python):/Python --volume $(abspath Bash):/Bash --tmpfs /tmpfs
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)"

# General
all:
.PHONY: all

log Output Docker/Tools:
	mkdir $@

# Tools
Docker/Tools/Mothur.zip: | Docker/Tools
	wget "https://github.com/mothur/mothur/releases/download/v1.44.2/Mothur.Ubuntu_18.zip" -O $@
TOOLS += Docker/Tools/Mothur.zip

# Docker
TOOLS += $(wildcard Docker/*)
build.log: Docker/Dockerfile $(TOOLS) | log Output
	rm -fv $@
	docker images | grep $(IMAGE_NAME) && docker rmi $(IMAGE_NAME) || true
	docker build --rm --tag $(IMAGE_NAME) $(<D) | tee $@

build: build.log
.PHONY: build

interactive: build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --interactive $(IMAGE_NAME) /bin/bash || true
.PHONY: interactive

delete: build.log
	docker rmi $(IMAGE_NAME)
	rm -fv build.log
.PHONY: delete

stop:
	docker rm $(CONTAINER_NAME)

# SGE
run: | log Output
	echo "make -C $(PWD) latest" > tmp.sh
	qsub -cwd -l h_vmem=$(MEMS) -m abe -M "230@fumire.moe" -N Premature_$(DATE) -pe smp $(CPUS) -o $(abspath log) -e $(abspath log) tmp.sh
.PHONY: run

# Actual
latest: step01 step02 step03 step04 step05 step06 step09 step21 step22
.PHONY: latest

# Step 01 (Make manifest files)
Output/Manifest:
	mkdir -p $@

Output/Manifest/first.tsv: Python/step01.py $(wildcard Data/RawData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/second.tsv: Python/step01.py $(wildcard Data/HN00145920/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/EBI.tsv: Python/step01-2.py $(wildcard Data/EBIData/ena_files/*/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/HMP.tsv: Python/step01-3.py $(wildcard Data/HMPData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/Stool.tsv: Python/step01-4.py $(wildcard Data/StoolData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

step01: Output/Manifest/first.tsv Output/Manifest/second.tsv Output/Manifest/EBI.tsv Output/Manifest/HMP.tsv Output/Manifest/Stool.tsv
.PHONY: step01

# Step 02 (De-multiplexing)
Output/Demux:
	mkdir -p $@

Output/Demux/first.qza: Output/Manifest/first.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/second.qza: Output/Manifest/second.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/EBI.qza: Output/Manifest/EBI.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[SequencesWithQuality]" --input-format "SingleEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/HMP.qza: Output/Manifest/HMP.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/Stool.qza: Output/Manifest/Stool.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/%.qzv: Output/Demux/%.qza | build.log
	$(DOCKER) qiime demux summarize --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step02: Output/Demux/first.qza Output/Demux/first.qzv Output/Demux/second.qza Output/Demux/second.qzv Output/Demux/EBI.qza Output/Demux/EBI.qzv Output/Demux/HMP.qza Output/Demux/HMP.qzv Output/Demux/Stool.qza Output/Demux/Stool.qzv
.PHONY: step02

# Step 03 (Quality filter)
Output/Filter:
	mkdir -p $@

Output/Filter/%.seq.qza Output/Filter/%.stat.qza: Output/Demux/%.qza | Output/Filter build.log
	$(DOCKER) qiime quality-filter q-score --i-demux $(addprefix /,$<) --o-filtered-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-filter-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step03: Output/Filter/first.seq.qza Output/Filter/second.seq.qza Output/Filter/EBI.seq.qza Output/Filter/HMP.seq.qza Output/Filter/Stool.seq.qza
.PHONY: step03

# Step 04 (DADA2 Denoising)
Output/DADA2:
	mkdir -p $@

Output/DADA2/first.table.qza Output/DADA2/first.seq.qza Output/DADA2/first.stat.qza: Output/Demux/first.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f 300 --p-trunc-len-r 265 --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/second.table.qza Output/DADA2/second.seq.qza Output/DADA2/second.stat.qza: Output/Demux/second.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f 300 --p-trunc-len-r 222 --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/EBI.table.qza Output/DADA2/EBI.seq.qza Output/DADA2/EBI.stat.qza: Output/Demux/EBI.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-single --i-demultiplexed-seqs $(addprefix /,$<) --p-trunc-len 150 --p-n-threads $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/HMP.table.qza Output/DADA2/HMP.seq.qza Output/DADA2/HMP.stat.qza: Output/Demux/HMP.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f 278 --p-trunc-len-r 226 --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/Stool.table.qza Output/DADA2/Stool.seq.qza Output/DADA2/Stool.stat.qza: Output/Demux/Stool.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f 250 --p-trunc-len-r 251 --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step04: Output/DADA2/first.seq.qza Output/DADA2/second.seq.qza Output/DADA2/EBI.seq.qza Output/DADA2/HMP.seq.qza Output/DADA2/Stool.seq.qza
.PHONY: step04

# Step 05 (Deblur Denoising)
Output/Deblur:
	mkdir -p $@

Output/Deblur/first.table.qza Output/Deblur/first.seq.qza Output/Deblur/first.stat.qza: Output/Filter/first.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length 265 --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/second.table.qza Output/Deblur/second.seq.qza Output/Deblur/second.stat.qza: Output/Filter/second.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length 222 --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/EBI.table.qza Output/Deblur/EBI.seq.qza Output/Deblur/EBI.stat.qza: Output/Filter/EBI.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length 150 --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/HMP.table.qza Output/Deblur/HMP.seq.qza Output/Deblur/HMP.stat.qza: Output/Filter/HMP.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length 226 --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/Stool.table.qza Output/Deblur/Stool.seq.qza Output/Deblur/Stool.stat.qza: Output/Filter/Stool.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length 250 --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step05: Output/Deblur/first.seq.qza Output/Deblur/second.seq.qza Output/Deblur/EBI.seq.qza Output/Deblur/HMP.seq.qza Output/Deblur/Stool.seq.qza
.PHONY: step05

# Step 06 (Taxonomy information)
Output/Taxonomy:
	mkdir -p $@

Output/Taxonomy/silva.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2020.6/common/silva-138-99-nb-classifier.qza" -O $@

Output/Taxonomy/gg.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2020.6/common/gg-13-8-99-nb-classifier.qza" -O $@

Output/Taxonomy/%.DADA2.silva.qza: Output/Taxonomy/silva.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.gg.qza: Output/Taxonomy/gg.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.silva.qza: Output/Taxonomy/silva.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.gg.qza: Output/Taxonomy/gg.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step06: Output/Taxonomy/first.DADA2.silva.qza Output/Taxonomy/first.DADA2.gg.qza Output/Taxonomy/first.Deblur.silva.qza Output/Taxonomy/first.Deblur.gg.qza Output/Taxonomy/second.DADA2.silva.qza Output/Taxonomy/second.DADA2.gg.qza Output/Taxonomy/second.Deblur.gg.qza Output/Taxonomy/second.Deblur.silva.qza Output/Taxonomy/EBI.DADA2.silva.qza Output/Taxonomy/EBI.DADA2.gg.qza Output/Taxonomy/EBI.Deblur.silva.qza Output/Taxonomy/EBI.Deblur.gg.qza Output/Taxonomy/HMP.DADA2.silva.qza Output/Taxonomy/HMP.DADA2.gg.qza Output/Taxonomy/HMP.Deblur.silva.qza Output/Taxonomy/HMP.Deblur.gg.qza
.PHONY: step06

# Step 07 (Build Metadata)
Output/Metadata:
	mkdir -p $@

Output/Metadata/first.tsv: Python/step07.py Data/Metadata/샘플정보.xlsx $(wildcard Data/firstData/*.fastq.gz) | Output/Metadata build.log
	@$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Metadata/EBI.tsv: Python/step07-2.py Data/EBIData/filereport_read_run_PRJEB10914.tsv $(wildcard Data/EBIData/ena_files/*/*.fastq.gz) | Output/Metadata build.log
	@$(DOCKER) python3 $(addprefix /,$^) > $@

step07: Output/Metadata/first.tsv Output/Metadata/EBI.tsv
.PHONY: step07

# Step 08 (Run ANCOM)
Output/ANCOM:
	mkdir -p $@

Output/ANCOM/%.DADA2.gg.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.DADA2.silva.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.gg.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.silva.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.CV.qza: Output/ANCOM/%.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[site] IN ('Cervix', 'Vagina')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.mouth.qza: Output/ANCOM/%.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[site] IN ('Mouth')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.pseudocount.qza: Output/ANCOM/%.qza | build.log
	$(DOCKER) qiime composition add-pseudocount --i-table $(addprefix /,$(word 1,$^)) --o-composition-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/first.%.site.qzv: Output/ANCOM/first.%.pseudocount.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "site" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/EBI.%.site.qzv: Output/ANCOM/EBI.%.pseudocount.qza Output/Metadata/EBI.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "site" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.site_premature.qzv: Output/ANCOM/%.pseudocount.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "site_premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.premature.qzv: Output/ANCOM/%.pseudocount.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "premature2" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.CV.premature.qzv: Output/ANCOM/%.CV.pseudocount.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "premature2" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.mouth.premature.qzv: Output/ANCOM/%.mouth.pseudocount.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "premature2" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%/ancom.tsv: Output/ANCOM/%.qzv | build.log
	$(DOCKER) rm -rfv $(addprefix /,$(@D))
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D))

step08: Output/ANCOM/first.DADA2.gg.site/ancom.tsv Output/ANCOM/first.DADA2.silva.site/ancom.tsv Output/ANCOM/first.Deblur.gg.site/ancom.tsv Output/ANCOM/first.Deblur.silva.site/ancom.tsv Output/ANCOM/EBI.DADA2.gg.site/ancom.tsv Output/ANCOM/EBI.DADA2.silva.site/ancom.tsv Output/ANCOM/EBI.Deblur.gg.site/ancom.tsv Output/ANCOM/EBI.Deblur.silva.site/ancom.tsv
.PHONY: step08

# Step 09 (Get TSV)
Output/TSV:
	mkdir -p $@

Output/TSV/first.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/first.table.qza Output/Taxonomy/first.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/first.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/first.table.qza Output/Taxonomy/first.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/second.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/second.table.qza Output/Taxonomy/second.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/second.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/second.table.qza Output/Taxonomy/second.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/Stool.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/Stool.table.qza Output/Taxonomy/Stool.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/Stool.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/Stool.table.qza Output/Taxonomy/Stool.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/EBI.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/EBI.table.qza Output/Taxonomy/EBI.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/EBI.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/EBI.table.qza Output/Taxonomy/EBI.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/HMP.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/HMP.table.qza Output/Taxonomy/HMP.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/HMP.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/HMP.table.qza Output/Taxonomy/HMP.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step09: Output/TSV/first.DADA2.gg.tsv Output/TSV/first.DADA2.silva.tsv Output/TSV/first.Deblur.gg.tsv Output/TSV/first.Deblur.silva.tsv Output/TSV/second.DADA2.gg.tsv  Output/TSV/second.DADA2.silva.tsv Output/TSV/second.Deblur.gg.tsv Output/TSV/second.Deblur.silva.tsv Output/TSV/Stool.DADA2.gg.tsv Output/TSV/Stool.DADA2.silva.tsv Output/TSV/Stool.Deblur.gg.tsv Output/TSV/Stool.Deblur.silva.tsv Output/TSV/EBI.DADA2.gg.tsv Output/TSV/EBI.DADA2.silva.tsv Output/TSV/EBI.Deblur.gg.tsv Output/TSV/EBI.Deblur.silva.tsv Output/TSV/HMP.DADA2.gg.tsv Output/TSV/HMP.DADA2.silva.tsv Output/TSV/HMP.Deblur.gg.tsv Output/TSV/HMP.Deblur.silva.tsv
.PHONY: step09

# Step 10 (Select & Clearify TSV)
Output/TSV/site.%.tar.gz: Python/step10.py Output/TSV/first.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^) $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/TSV/premature.%.tar.gz: Python/step10.py Output/TSV/first.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^) $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step10: 
.PHONY: step10

# Step 11 (Making t-SNE)
Output/TSNE:
	mkdir -p $@

Output/TSNE/%.tar.gz: Python/step11.py Output/TSV/%.tar.gz | Output/TSNE build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

step11: Output/TSNE/site.DADA2.gg.tar.gz Output/TSNE/site.DADA2.silva.tar.gz Output/TSNE/site.Deblur.gg.tar.gz Output/TSNE/site.Deblur.silva.tar.gz
.PHONY: step11

# Step 12 (Draw t-SNE)
Output/TSNE/%.png: Python/step12.py Output/TSNE/%.tar.gz Output/Metadata/first.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step12: Output/TSNE/site.DADA2.gg.png Output/TSNE/premature.DADA2.gg.png Output/TSNE/site.DADA2.silva.png Output/TSNE/premature.DADA2.silva.png Output/TSNE/site.Deblur.gg.png Output/TSNE/premature.Deblur.gg.png Output/TSNE/site.Deblur.silva.png Output/TSNE/premature.Deblur.silva.png
.PHONY: step12

# Step 21 (Read & clearify first TSV into pandas)
Output/TSV/first.%.tar.gz: Python/step21.py Output/TSV/first.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/second.%.tar.gz: Python/step21.py Output/TSV/second.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/Stool.%.tar.gz: Python/step21.py Output/TSV/Stool.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/EBI.%.tar.gz: Python/step21.py Output/TSV/EBI.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/HMP.%.tar.gz: Python/step21.py Output/TSV/HMP.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step21: Output/TSV/first.DADA2.gg.tar.gz Output/TSV/first.DADA2.silva.tar.gz Output/TSV/first.Deblur.gg.tar.gz Output/TSV/first.Deblur.silva.tar.gz Output/TSV/second.DADA2.gg.tar.gz Output/TSV/second.DADA2.silva.tar.gz Output/TSV/second.Deblur.gg.tar.gz Output/TSV/second.Deblur.silva.tar.gz Output/TSV/Stool.DADA2.gg.tar.gz Output/TSV/Stool.DADA2.silva.tar.gz Output/TSV/Stool.Deblur.gg.tar.gz Output/TSV/Stool.Deblur.silva.tar.gz Output/TSV/HMP.DADA2.gg.tar.gz Output/TSV/HMP.DADA2.silva.tar.gz Output/TSV/HMP.Deblur.gg.tar.gz Output/TSV/HMP.Deblur.silva.tar.gz Output/TSV/EBI.DADA2.gg.tar.gz Output/TSV/EBI.DADA2.silva.tar.gz Output/TSV/EBI.Deblur.gg.tar.gz Output/TSV/EBI.Deblur.silva.tar.gz
.PHONY: step21

# Step 22 (Making & Draw t-SNE)
Output/Step22:
	mkdir -p $@

Output/Step22/%.tar.gz: Python/step11.py Output/TSV/%.tar.gz | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step22/%.png: Python/step22.py Output/Step22/%.tar.gz | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step22/DB.%.pdf: Python/step22-2.py Output/TSV/first.%.tar.gz Output/TSV/second.%.tar.gz Output/TSV/Stool.%.tar.gz Output/TSV/EBI.%.tar.gz Output/TSV/HMP.%.tar.gz | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

step22: Output/Step22/DB.DADA2.gg.pdf Output/Step22/DB.DADA2.silva.pdf Output/Step22/DB.Deblur.gg.pdf Output/Step22/DB.Deblur.silva.pdf
.PHONY: step22

# Step 23 (Draw histogram)
Output/Step23:
	mkdir -p $@

Output/Step23/age.png: Python/step23.py Data/Metadata/샘플정보.xlsx | Output/Step23 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "나이"

Output/Step23/weight.png: Python/step23.py Data/Metadata/샘플정보.xlsx | Output/Step23 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "임신전 몸무게"

Output/Step23/BMI.png: Python/step23.py Data/Metadata/샘플정보.xlsx | Output/Step23 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "BMI"

step23: Output/Step23/age.png Output/Step23/weight.png Output/Step23/BMI.png
.PHONY: step23

# Step 24 (Draw histogram with parsing)
Output/Step24:
	mkdir -p $@

Output/Step24/weight_now.png: Python/step24.py Data/Metadata/샘플정보.xlsx | Output/Step24 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "체중/키" 0 "/" 1> $@.stdout 2> $@.stderr

Output/Step24/height.png: Python/step24.py Data/Metadata/샘플정보.xlsx | Output/Step24 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "체중/키" 1 "/" 1> $@.stdout 2> $@.stderr

Output/Step24/weeks.png: Python/step24.py Data/Metadata/샘플정보.xlsx | Output/Step24 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "분만주수" 0 "+" 1> $@.stdout 2> $@.stderr

step24: Output/Step24/weight_now.png Output/Step24/height.png Output/Step24/weeks.png
.PHONY: step24

# Step 26 (RandomForest Classifier)
Output/Step26:
	mkdir -p $@

Output/Step26/RF.%.tar: Python/step26.py Output/TSV/helixco.%.tar.gz Output/TSV/EBI.%.tar.gz Output/TSV/HMP.%.tar.gz Output/Metadata/first.tsv | Output/Step26 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

step26: Output/Step26/RF.DADA2.gg.tar Output/Step26/RF.DADA2.silva.tar Output/Step26/RF.Deblur.gg.tar Output/Step26/RF.Deblur.silva.tar
.PHONY: step26

# Step 27 (DecisionTree Classifier)
Output/Step27:
	mkdir -p $@

Output/Step27/DT.%.tar: Python/step27.py Output/TSV/helixco.%.tar.gz Output/TSV/EBI.%.tar.gz Output/TSV/HMP.%.tar.gz Output/Metadata/first.tsv | Output/Step27 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

step27: Output/Step27/DT.DADA2.gg.tar Output/Step27/DT.DADA2.silva.tar Output/Step27/DT.Deblur.gg.tar Output/Step27/DT.Deblur.silva.tar
.PHONY: step27
