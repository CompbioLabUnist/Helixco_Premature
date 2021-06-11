# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = premature:latest
CPUS = 20
MEMS = 300G
PWD := $(shell pwd)
TOOLS = 
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python --volume $(abspath Bash):/Bash:ro --tmpfs /tmpfs
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
latest: step01 step02 step03 step04 step05 step06# step07 step08 step09 step21 step22 step26 step28 step29 step30
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

Output/Manifest/Stool.tsv: Python/step01-4.py Data/Metadata/readable.xlsx $(wildcard Data/StoolData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$(filter %.py,$^) $(filter %.xlsx,$^)) $(foreach f,$(wildcard Data/StoolData/*.fastq.gz),--stool $(addprefix /,$f)) > $@

Output/Manifest/everything.tsv: Python/step01-5.py Data/Metadata/readable.xlsx $(wildcard Data/RawData/*.fastq.gz) $(wildcard Data/HN00145920/*.fastq.gz) $(wildcard Data/HN00148220/*.fastq.gz) $(wildcard Data/StoolData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$(filter %.py,$^) $(filter %.xlsx,$^)) $(foreach f,$(wildcard Data/RawData/*.fastq.gz),--first $(addprefix /,$f)) $(foreach f,$(wildcard Data/HN00145920/*.fastq.gz),--second $(addprefix /,$f)) $(foreach f,$(wildcard Data/HN00148220/*.fastq.gz),--third $(addprefix /,$f)) > $@

step01: Output/Manifest/Stool.tsv Output/Manifest/everything.tsv
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

Output/Demux/everything.qza: Output/Manifest/everything.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/%.qzv: Output/Demux/%.qza | build.log
	$(DOCKER) qiime demux summarize --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step02: Output/Demux/Stool.qza Output/Demux/Stool.qzv Output/Demux/everything.qza Output/Demux/everything.qzv
.PHONY: step02

# Step 03 (Quality filter)
Output/Filter:
	mkdir -p $@

Output/Filter/%.seq.qza Output/Filter/%.stat.qza: Output/Demux/%.qza | Output/Filter build.log
	$(DOCKER) qiime quality-filter q-score --i-demux $(addprefix /,$<) --o-filtered-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-filter-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step03: Output/Filter/Stool.seq.qza Output/Filter/everything.seq.qza
.PHONY: step03

# Step 04 (DADA2 Denoising)
Output/DADA2:
	mkdir -p $@

Output/DADA2/first.table.qza Output/DADA2/first.seq.qza Output/DADA2/first.stat.qza: Output/Demux/first.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "265" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/second.table.qza Output/DADA2/second.seq.qza Output/DADA2/second.stat.qza: Output/Demux/second.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "222" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/EBI.table.qza Output/DADA2/EBI.seq.qza Output/DADA2/EBI.stat.qza: Output/Demux/EBI.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-single --i-demultiplexed-seqs $(addprefix /,$<) --p-trunc-len "150" --p-n-threads $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/HMP.table.qza Output/DADA2/HMP.seq.qza Output/DADA2/HMP.stat.qza: Output/Demux/HMP.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "278" --p-trunc-len-r "226" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/Stool.table.qza Output/DADA2/Stool.seq.qza Output/DADA2/Stool.stat.qza: Output/Demux/Stool.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "250" --p-trunc-len-r "251" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/everything.table.qza Output/DADA2/everything.seq.qza Output/DADA2/everything.stat.qza: Output/Demux/everything.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "245" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step04: Output/DADA2/Stool.seq.qza Output/DADA2/everything.seq.qza
.PHONY: step04

# Step 05 (Deblur Denoising)
Output/Deblur:
	mkdir -p $@

Output/Deblur/first.table.qza Output/Deblur/first.seq.qza Output/Deblur/first.stat.qza: Output/Filter/first.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "265" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/second.table.qza Output/Deblur/second.seq.qza Output/Deblur/second.stat.qza: Output/Filter/second.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "222" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/EBI.table.qza Output/Deblur/EBI.seq.qza Output/Deblur/EBI.stat.qza: Output/Filter/EBI.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "150" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/HMP.table.qza Output/Deblur/HMP.seq.qza Output/Deblur/HMP.stat.qza: Output/Filter/HMP.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "226" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/Stool.table.qza Output/Deblur/Stool.seq.qza Output/Deblur/Stool.stat.qza: Output/Filter/Stool.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "250" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/everything.table.qza Output/Deblur/everything.seq.qza Output/Deblur/everything.stat.qza: Output/Filter/everything.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "245" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step05: Output/Deblur/Stool.seq.qza Output/Deblur/everything.seq.qza
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

step06: Output/Taxonomy/everything.DADA2.gg.qza Output/Taxonomy/everything.DADA2.silva.qza Output/Taxonomy/everything.Deblur.gg.qza Output/Taxonomy/everything.Deblur.silva.qza Output/Taxonomy/Stool.DADA2.gg.qza Output/Taxonomy/Stool.DADA2.silva.qza Output/Taxonomy/Stool.Deblur.gg.qza Output/Taxonomy/Stool.Deblur.silva.qza
.PHONY: step06

# Step 07 (Build Metadata)
Output/Metadata:
	mkdir -p $@

Output/Metadata/first.tsv: Python/step07.py Data/Metadata/sample_information.xlsx $(sort $(wildcard Data/RawData/*.fastq.gz)) | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Metadata/EBI.tsv: Python/step07-2.py Data/EBIData/filereport_read_run_PRJEB10914.tsv $(sort $(wildcard Data/EBIData/ena_files/*/*.fastq.gz)) | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^)

Output/Metadata/second.tsv: Python/step07-3.py Data/Metadata/sample_information.xlsx $(sort $(wildcard Data/HN00145920/*.fastq.gz)) | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^)

Output/Metadata/everything.tsv: Python/step07-4.py Output/Manifest/everything.tsv Data/Metadata/metadata.xlsx | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Metadata/Stool.tsv: Python/step07-4.py Output/Manifest/Stoo.tsv Data/Metadata/metadata.xlsx | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

step07: Output/Metadata/Stool.tsv Output/Metadata/everything.tsv
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

Output/ANCOM/everything.%.Mouth.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Mouth'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Cervix.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Cervix'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Vagina.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Vagina'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Baby.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site] IN ('Baby-1day', 'Baby-3day', 'Baby-5day')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Stool.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site] IN ('Stool-1day', 'Stool-3day', 'Stool-5day')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.pseudocount.qza: Output/ANCOM/%.qza | build.log
	$(DOCKER) qiime composition add-pseudocount --i-table $(addprefix /,$(word 1,$^)) --o-composition-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/first.%.Site.qzv: Output/ANCOM/first.%.pseudocount.qza Output/Metadata/first.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Site" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.PTB.qzv: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Detail_PTB.qzv: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Detail Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%/ancom.tsv Output/ANCOM/%/data.tsv: Output/ANCOM/%.qzv | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D))
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

Output/ANCOM/%.pdf: Python/step08.py Output/ANCOM/%/ancom.tsv Output/ANCOM/%/data.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step08-DADA2-gg: Output/ANCOM/everything.DADA2.gg.qza Output/ANCOM/everything.DADA2.gg.Mouth.qza Output/ANCOM/everything.DADA2.gg.Cervix.qza Output/ANCOM/everything.DADA2.gg.Vagina.qza Output/ANCOM/everything.DADA2.gg.Baby.qza Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.Detail_PTB.pdf
step08-DADA2-silva: 
step08-Deblur-gg: 
step08-Deblur-silva: 
step08: step08-DADA2-gg step08-DADA2-silva step08-Deblur-gg step08-Deblur-silva
.PHONY: step08 step08-DADA2-gg step08-DADA2-silva step08-Deblur-gg step08-Deblur-silva

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

Output/TSV/everything.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/everything.table.qza Output/Taxonomy/everything.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/everything.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/everything.table.qza Output/Taxonomy/everything.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step09: Output/TSV/everything.DADA2.gg.tsv Output/TSV/everything.DADA2.silva.tsv Output/TSV/everything.Deblur.gg.tsv Output/TSV/everything.Deblur.silva.tsv
.PHONY: step09

# Step 21 (Read & clearify first TSV into pandas)
Output/TSV/everything.%.merge.tsv: Python/step21.py Output/TSV/everything.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step21: Output/TSV/everything.DADA2.gg.merge.tsv Output/TSV/everything.DADA2.silva.merge.tsv Output/TSV/everything.Deblur.gg.merge.tsv Output/TSV/everything.Deblur.silva.merge.tsv
.PHONY: step21

# Step 22 (Making & Draw t-SNE)
Output/Step22:
	mkdir -p $@

Output/Step22/everything.%.tar.gz: Python/step11.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step22/%.tar: Python/step22-3.py Output/Step22/%.tar.gz | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step22: Output/Step22/everything.DADA2.gg.tar.gz Output/Step22/everything.DADA2.gg.tar Output/Step22/everything.DADA2.silva.tar.gz Output/Step22/everything.DADA2.silva.tar Output/Step22/everything.Deblur.gg.tar.gz Output/Step22/everything.Deblur.gg.tar Output/Step22/everything.Deblur.silva.tar.gz Output/Step22/everything.Deblur.silva.tar
.PHONY: step22

# Step 26 (RandomForest Classifier)
Output/Step26:
	mkdir -p $@

Output/Step26/RF.%.tar: Python/step26.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | Output/Step26 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step26: Output/Step26/RF.DADA2.gg.tar Output/Step26/RF.DADA2.silva.tar Output/Step26/RF.Deblur.gg.tar Output/Step26/RF.Deblur.silva.tar
.PHONY: step26

# Step 27 (DecisionTree Classifier)
Output/Step27:
	mkdir -p $@

Output/Step27/DT.%.tar: Python/step27.py Output/TSV/helixco.%.tar.gz Output/TSV/EBI.%.tar.gz Output/TSV/HMP.%.tar.gz Output/Metadata/first.tsv | Output/Step27 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

step27: Output/Step27/DT.DADA2.gg.tar Output/Step27/DT.DADA2.silva.tar Output/Step27/DT.Deblur.gg.tar Output/Step27/DT.Deblur.silva.tar
.PHONY: step27

# Step 28 (Phylogenetic diversity analysis)
Output/Step28:
	mkdir -p $@

Output/Step28/%.DADA2.aligned-rep-seqs.qza Output/Step28/%.DADA2.masked-aligned-rep-seqs.qza Output/Step28/%.DADA2.unrooted-tree.qza Output/Step28/%.DADA2.rooted-tree.qza: Output/DADA2/%.seq.qza | Output/Step28 build.log
	$(DOCKER) qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $(addprefix /,$<) --o-alignment $(addprefix /,$(basename $(basename $@))).aligned-rep-seqs.qza --o-masked-alignment $(addprefix /,$(basename $(basename $@))).masked-aligned-rep-seqs.qza --o-tree $(addprefix /,$(basename $(basename $@))).unrooted-tree.qza --o-rooted-tree $(addprefix /,$(basename $(basename $@))).rooted-tree.qza --p-n-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step28/%.Deblur.aligned-rep-seqs.qza Output/Step28/%.Deblur.masked-aligned-rep-seqs.qza Output/Step28/%.Deblur.unrooted-tree.qza Output/Step28/%.Deblur.rooted-tree.qza: Output/Deblur/%.seq.qza | Output/Step28 build.log
	$(DOCKER) qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $(addprefix /,$<) --o-alignment $(addprefix /,$(basename $(basename $@))).aligned-rep-seqs.qza --o-masked-alignment $(addprefix /,$(basename $(basename $@))).masked-aligned-rep-seqs.qza --o-tree $(addprefix /,$(basename $(basename $@))).unrooted-tree.qza --o-rooted-tree $(addprefix /,$(basename $(basename $@))).rooted-tree.qza --p-n-threads $(CPUS) 1> $@.stdout 2> $@.stderr

step28: Output/Step28/everything.DADA2.rooted-tree.qza
.PHONY: step28

# Step 29 (Core phylogenetic metrics)
Output/Step29 Output/Step29/DADA2 Output/Step29/Deblur:
	mkdir -p $@

Output/Step29/%.DADA2.table.qzv: Output/DADA2/%.table.qza Output/Metadata/%.tsv | Output/Step29 build.log
	$(DOCKER) qiime feature-table summarize --i-table $(addprefix /,$(word 1,$^)) --m-sample-metadata-file $(addprefix /,$(word 2,$^)) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.DADA2.seq.qzv: Output/DADA2/%.seq.qza | Output/Step29 build.log
	$(DOCKER) qiime feature-table tabulate-seqs --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.Deblur.table.qzv: Output/Deblur/%.table.qza Output/Metadata/%.tsv | Output/Step29 build.log
	$(DOCKER) qiime feature-table summarize --i-table $(addprefix /,$(word 1,$^)) --m-sample-metadata-file $(addprefix /,$(word 2,$^)) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.Deblur.seq.qzv: Output/Deblur/%.seq.qza | Output/Step29 build.log
	$(DOCKER) qiime feature-table tabulate-seqs --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/DADA2/everything.faith_pd_vector.qza Output/Step29/DADA2/everything.unweighted_unifrac_distance_matrix.qza Output/Step29/DADA2/everything.bray_curtis_pcoa_results.qza Output/Step29/DADA2/everything.shannon_vector.qza Output/Step29/DADA2/everything.rarefied_table.qza Output/Step29/DADA2/everything.weighted_unifrac_distance_matrix.qza Output/Step29/DADA2/everything.jaccard_pcoa_results.qza Output/Step29/DADA2/everything.weighted_unifrac_pcoa_results.qza Output/Step29/DADA2/everything.observed_features_vector.qza Output/Step29/DADA2/everything.jaccard_distance_matrix.qza Output/Step29/DADA2/everything.evenness_vector.qza Output/Step29/DADA2/everything.bray_curtis_distance_matrix.qza Output/Step29/DADA2/everything.unweighted_unifrac_pcoa_results.qza Output/Step29/DADA2/everything.unweighted_unifrac_emperor.qzv Output/Step29/DADA2/everything.jaccard_emperor.qzv Output/Step29/DADA2/everything.bray_curtis_emperor.qzv Output/Step29/DADA2/everything.weighted_unifrac_emperor.qzv: Output/Step28/everything.DADA2.rooted-tree.qza Output/DADA2/everything.table.qza Output/Metadata/everything.tsv | Output/Step29/DADA2 build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "1046" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step29/Deblur/everything.faith_pd_vector.qza Output/Step29/Deblur/everything.unweighted_unifrac_distance_matrix.qza Output/Step29/Deblur/everything.bray_curtis_pcoa_results.qza Output/Step29/Deblur/everything.shannon_vector.qza Output/Step29/Deblur/everything.rarefied_table.qza Output/Step29/Deblur/everything.weighted_unifrac_distance_matrix.qza Output/Step29/Deblur/everything.jaccard_pcoa_results.qza Output/Step29/Deblur/everything.weighted_unifrac_pcoa_results.qza Output/Step29/Deblur/everything.observed_features_vector.qza Output/Step29/Deblur/.jaccard_distance_matrix.qza Output/Step29/Deblur/everything.evenness_vector.qza Output/Step29/Deblur/everything.bray_curtis_distance_matrix.qza Output/Step29/Deblur/everything.unweighted_unifrac_pcoa_results.qza Output/Step29/Deblur/everything.unweighted_unifrac_emperor.qzv Output/Step29/Deblur/everything.jaccard_emperor.qzv Output/Step29/Deblur/everything.bray_curtis_emperor.qzv Output/Step29/Deblur/everything.weighted_unifrac_emperor.qzv: Output/Step28/everything.Deblur.rooted-tree.qza Output/Deblur/everything.table.qza Output/Metadata/everything.tsv | Output/Step29/Deblur build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "4864" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

step29-depth: Output/Step29/everything.DADA2.table.qzv Output/Step29/everything.DADA2.seq.qzv Output/Step29/everything.Deblur.table.qzv Output/Step29/everything.Deblur.seq.qzv
step29: Output/Step29/DADA2/everything.faith_pd_vector.qza Output/Step29/Deblur/everything.faith_pd_vector.qza
.PHONY: step29-depth step29

# Step 30 (Alpha-diversity)
Output/Step30/DADA2 Output/Step30/Deblur:
	mkdir -p $@

Output/Step30/DADA2/%.tsv: Output/Step29/DADA2/%_vector.qza | Output/Step30/DADA2 build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/alpha-diversity.tsv $(addprefix /,$@)

Output/Step30/Deblur/%.tsv: Output/Step29/Deblur/%_vector.qza | Output/Step30/Deblur build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/alpha-diversity.tsv $(addprefix /,$@)

step30: Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv
.PHONY: step30

# Step 31 (Violin plots from ANCOM)
Output/ANCOM/everything.DADA2.gg.%.tar: Python/step31.py Output/ANCOM/everything.DADA2.gg.%/ancom.tsv Output/TSV/everything.DADA2.gg.merge.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@)

step31: Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.tar# Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.Detail_PTB.tar
